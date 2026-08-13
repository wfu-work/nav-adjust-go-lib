package sparse

import (
	"context"
	"errors"
	"fmt"
	"math"
)

var (
	// ErrBreakdown indicates that the matrix is not numerically positive
	// definite for the conjugate-gradient iteration.
	ErrBreakdown = errors.New("sparse: conjugate-gradient breakdown")
	// ErrNotConverged indicates that the configured iteration limit was reached.
	ErrNotConverged = errors.New("sparse: conjugate-gradient did not converge")
)

// PCGOptions controls preconditioned conjugate-gradient iteration.
type PCGOptions struct {
	MaxIterations     int
	RelativeTolerance float64
	AbsoluteTolerance float64
	Preconditioner    Preconditioner
}

// PCGInfo describes the achieved convergence.
type PCGInfo struct {
	Iterations       int
	ResidualNorm     float64
	RelativeResidual float64
}

// ProjectFunc projects a vector into the linear subspace in which PCG should
// iterate. Implementations may reuse private scratch storage because calls are
// synchronous within one solve.
type ProjectFunc func(vector []float64) error

// SolvePCG solves matrix*x=rhs with the configured preconditioner, defaulting
// to scalar Jacobi. The initial solution is zero and the returned slice never
// aliases rhs.
func SolvePCG(matrix *Matrix, rhs []float64, options PCGOptions) ([]float64, PCGInfo, error) {
	return SolvePCGContext(context.Background(), matrix, rhs, options)
}

// SolvePCGContext is SolvePCG with cancellation support.
func SolvePCGContext(ctx context.Context, matrix *Matrix, rhs []float64, options PCGOptions) ([]float64, PCGInfo, error) {
	return solvePCG(ctx, matrix, rhs, nil, options)
}

// SolveProjectedPCG solves the symmetric system in the subspace selected by
// project. Both the right-hand side and every preconditioned search direction
// are projected, so the result remains in that subspace. This is used for
// equality-constrained adjustment without forming an indefinite KKT matrix.
func SolveProjectedPCG(matrix *Matrix, rhs []float64, project ProjectFunc, options PCGOptions) ([]float64, PCGInfo, error) {
	return SolveProjectedPCGContext(context.Background(), matrix, rhs, project, options)
}

// SolveProjectedPCGContext is SolveProjectedPCG with cancellation support.
func SolveProjectedPCGContext(ctx context.Context, matrix *Matrix, rhs []float64, project ProjectFunc, options PCGOptions) ([]float64, PCGInfo, error) {
	if project == nil {
		return nil, PCGInfo{}, fmt.Errorf("sparse: projected PCG requires a projector")
	}
	return solvePCG(ctx, matrix, rhs, project, options)
}

func solvePCG(ctx context.Context, matrix *Matrix, rhs []float64, project ProjectFunc, options PCGOptions) ([]float64, PCGInfo, error) {
	if ctx == nil {
		return nil, PCGInfo{}, fmt.Errorf("sparse: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, PCGInfo{}, err
	}
	if matrix == nil || len(rhs) != matrix.Size() {
		return nil, PCGInfo{}, fmt.Errorf("sparse: PCG dimension mismatch")
	}
	n := matrix.Size()
	if options.MaxIterations <= 0 {
		options.MaxIterations = max(100, 10*n)
	}
	if options.RelativeTolerance <= 0 {
		options.RelativeTolerance = 1e-10
	}
	if options.AbsoluteTolerance <= 0 {
		options.AbsoluteTolerance = 1e-12
	}

	x := make([]float64, n)
	r := append([]float64(nil), rhs...)
	if project != nil {
		if err := project(r); err != nil {
			return nil, PCGInfo{}, fmt.Errorf("sparse: project right-hand side: %w", err)
		}
	}
	z := make([]float64, n)
	p := make([]float64, n)
	ap := make([]float64, n)
	preconditioner := options.Preconditioner
	if preconditioner == nil {
		var err error
		preconditioner, err = NewJacobiPreconditioner(matrix, project)
		if err != nil {
			return nil, PCGInfo{}, err
		}
	}
	rhsNorm := norm(r)
	tolerance := math.Max(options.AbsoluteTolerance, options.RelativeTolerance*rhsNorm)
	residualNorm := norm(r)
	info := PCGInfo{ResidualNorm: residualNorm}
	if rhsNorm > 0 {
		info.RelativeResidual = residualNorm / rhsNorm
	}
	if residualNorm <= tolerance {
		return x, info, nil
	}
	if err := preconditioner.Apply(z, r); err != nil {
		return nil, info, err
	}
	if project != nil {
		if err := project(z); err != nil {
			return nil, info, fmt.Errorf("sparse: project preconditioned residual: %w", err)
		}
	}
	copy(p, z)
	rz := dot(r, z)
	if rz <= 0 || !finite(rz) {
		return nil, info, fmt.Errorf("%w: invalid preconditioned residual", ErrBreakdown)
	}

	for iteration := 1; iteration <= options.MaxIterations; iteration++ {
		if err := ctx.Err(); err != nil {
			return nil, info, err
		}
		matrix.MulVec(ap, p)
		if project != nil {
			if err := project(ap); err != nil {
				return nil, info, fmt.Errorf("sparse: project matrix product: %w", err)
			}
		}
		curvature := dot(p, ap)
		if curvature <= 0 || !finite(curvature) {
			return nil, info, fmt.Errorf("%w: non-positive curvature at iteration %d", ErrBreakdown, iteration)
		}
		alpha := rz / curvature
		for i := range n {
			x[i] += alpha * p[i]
			r[i] -= alpha * ap[i]
		}
		if project != nil {
			if err := project(r); err != nil {
				return nil, info, fmt.Errorf("sparse: project residual: %w", err)
			}
		}
		residualNorm = norm(r)
		info.Iterations = iteration
		info.ResidualNorm = residualNorm
		if rhsNorm > 0 {
			info.RelativeResidual = residualNorm / rhsNorm
		}
		if !finite(residualNorm) {
			return nil, info, fmt.Errorf("%w: non-finite residual at iteration %d", ErrBreakdown, iteration)
		}
		if residualNorm <= tolerance {
			return x, info, nil
		}
		if err := preconditioner.Apply(z, r); err != nil {
			return nil, info, err
		}
		if project != nil {
			if err := project(z); err != nil {
				return nil, info, fmt.Errorf("sparse: project preconditioned residual: %w", err)
			}
		}
		nextRZ := dot(r, z)
		if nextRZ <= 0 || !finite(nextRZ) {
			return nil, info, fmt.Errorf("%w: invalid residual at iteration %d", ErrBreakdown, iteration)
		}
		beta := nextRZ / rz
		for i := range n {
			p[i] = z[i] + beta*p[i]
		}
		rz = nextRZ
	}
	return nil, info, fmt.Errorf("%w after %d iterations: relative residual %.6g", ErrNotConverged, info.Iterations, info.RelativeResidual)
}

func dot(left, right []float64) float64 {
	value := 0.0
	for i := range left {
		value += left[i] * right[i]
	}
	return value
}

func norm(values []float64) float64 { return math.Sqrt(math.Max(0, dot(values, values))) }

func finite(value float64) bool { return !math.IsNaN(value) && !math.IsInf(value, 0) }
