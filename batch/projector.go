package batch

import (
	"context"
	"errors"
	"fmt"
	"math"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
	"github.com/wfu-work/nav-adjust-go-lib/internal/sparse"
	"gonum.org/v1/gonum/mat"
)

// constraintProjector represents the Euclidean null-space projector
// P=I-C^T(CC^T)^-1C without materializing either C or P as a dense matrix.
// Only the small constraint Gram matrix is factorized densely.
type constraintProjector struct {
	rows       []exactRow
	parameters int
	factor     mat.Cholesky
}

func newConstraintProjector(compiled compiledProblem) (*constraintProjector, error) {
	count := len(compiled.exact)
	if count == 0 {
		return nil, invalid("projected solve requires at least one exact constraint")
	}
	gram := mat.NewSymDense(count, nil)
	for row := range count {
		for column := 0; column <= row; column++ {
			gram.SetSym(row, column, exactRowDot(compiled.exact[row], compiled.exact[column]))
		}
	}
	projector := &constraintProjector{rows: compiled.exact, parameters: compiled.parameters}
	if !projector.factor.Factorize(gram) {
		return nil, fmt.Errorf("%w: exact constraints are linearly dependent", adjust.ErrRankDeficient)
	}
	condition := projector.factor.Cond()
	if math.IsInf(condition, 1) || condition > 1e16 {
		return nil, fmt.Errorf("%w: exact-constraint Gram matrix has condition number %.6g", adjust.ErrRankDeficient, condition)
	}
	return projector, nil
}

func exactRowDot(left, right exactRow) float64 {
	value := 0.0
	leftIndex, rightIndex := 0, 0
	for leftIndex < len(left.terms) && rightIndex < len(right.terms) {
		leftTerm, rightTerm := left.terms[leftIndex], right.terms[rightIndex]
		switch {
		case leftTerm.Parameter < rightTerm.Parameter:
			leftIndex++
		case leftTerm.Parameter > rightTerm.Parameter:
			rightIndex++
		default:
			value += leftTerm.Coefficient * rightTerm.Coefficient
			leftIndex++
			rightIndex++
		}
	}
	return value
}

func (projector *constraintProjector) projectFunc() sparse.ProjectFunc {
	rhs := mat.NewVecDense(len(projector.rows), nil)
	multipliers := mat.NewVecDense(len(projector.rows), nil)
	return func(vector []float64) error {
		if len(vector) != projector.parameters {
			return fmt.Errorf("constraint projection dimension mismatch")
		}
		for row, constraint := range projector.rows {
			rhs.SetVec(row, dotTerms(constraint.terms, vector))
		}
		if err := projector.factor.SolveVecTo(multipliers, rhs); err != nil {
			return fmt.Errorf("solve exact-constraint Gram system: %w", err)
		}
		for row, constraint := range projector.rows {
			multiplier := multipliers.AtVec(row)
			for _, term := range constraint.terms {
				vector[term.Parameter] -= term.Coefficient * multiplier
			}
		}
		return nil
	}
}

func (projector *constraintProjector) particularSolution() ([]float64, error) {
	values := mat.NewVecDense(len(projector.rows), nil)
	for row, constraint := range projector.rows {
		values.SetVec(row, constraint.value)
	}
	multipliers := mat.NewVecDense(len(projector.rows), nil)
	if err := projector.factor.SolveVecTo(multipliers, values); err != nil {
		return nil, fmt.Errorf("batch: exact-constraint particular solution: %w", err)
	}
	result := make([]float64, projector.parameters)
	for row, constraint := range projector.rows {
		multiplier := multipliers.AtVec(row)
		for _, term := range constraint.terms {
			result[term.Parameter] += term.Coefficient * multiplier
		}
	}
	return result, nil
}

func solveProjectedSparse(ctx context.Context, normal sparseNormalEquation, compiled compiledProblem, batchOptions Options, options sparse.PCGOptions) (
	delta []float64,
	provider covarianceProvider,
	iterations int,
	relativeResidual float64,
	err error,
) {
	projector, err := newConstraintProjector(compiled)
	if err != nil {
		return nil, nil, 0, 0, err
	}
	options.Preconditioner, err = makePreconditioner(normal.n, batchOptions, projector.projectFunc())
	if err != nil {
		return nil, nil, 0, 0, mapPCGError("preconditioner", err)
	}
	particular, err := projector.particularSolution()
	if err != nil {
		return nil, nil, 0, 0, err
	}
	normalTimesParticular := make([]float64, compiled.parameters)
	normal.n.MulVec(normalTimesParticular, particular)
	rhs := make([]float64, compiled.parameters)
	for i := range rhs {
		rhs[i] = normal.u[i] - normalTimesParticular[i]
	}
	correction, info, err := sparse.SolveProjectedPCGContext(ctx, normal.n, rhs, projector.projectFunc(), options)
	if err != nil {
		return nil, nil, info.Iterations, info.RelativeResidual, mapPCGError("projected normal equation", err)
	}
	for i := range particular {
		particular[i] += correction[i]
	}
	provider = &projectedPCGCovarianceProvider{matrix: normal.n, options: options, projector: projector}
	return particular, provider, info.Iterations, info.RelativeResidual, nil
}

type projectedPCGCovarianceProvider struct {
	matrix    *sparse.Matrix
	options   sparse.PCGOptions
	projector *constraintProjector
}

func (provider *projectedPCGCovarianceProvider) column(ctx context.Context, index int) ([]float64, error) {
	rhs := make([]float64, provider.matrix.Size())
	rhs[index] = 1
	column, _, err := sparse.SolveProjectedPCGContext(ctx, provider.matrix, rhs, provider.projector.projectFunc(), provider.options)
	if err != nil {
		if errors.Is(err, sparse.ErrNotConverged) || errors.Is(err, sparse.ErrBreakdown) {
			return nil, mapPCGError("constrained covariance query", err)
		}
		return nil, fmt.Errorf("batch: constrained covariance query: %w", err)
	}
	return column, nil
}
