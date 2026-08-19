// Package batch implements linear weighted least-squares adjustment.
package batch

import "math"

// SolverMethod selects the normal-equation solver.
type SolverMethod string

// PreconditionerMethod selects the sparse PCG preconditioner.
type PreconditionerMethod string

const (
	// SolverDense uses dense Cholesky or constrained KKT factorization.
	SolverDense SolverMethod = "dense"
	// SolverAuto uses dense factorization for full covariance or small systems,
	// and sparse PCG for larger systems with reduced covariance output.
	SolverAuto SolverMethod = "auto"
	// SolverPCG uses a sparse normal matrix and preconditioned conjugate
	// gradients. Exact constraints are handled by null-space projection.
	SolverPCG SolverMethod = "pcg"
)

const (
	PreconditionerJacobi      PreconditionerMethod = "jacobi"
	PreconditionerBlockJacobi PreconditionerMethod = "block-jacobi"
	PreconditionerIC0         PreconditionerMethod = "ic0"
)

// CovarianceMode controls whether Solve materializes the complete inverse
// normal matrix. CovarianceNone still permits on-demand column queries through
// SolveDetailed.
type CovarianceMode string

const (
	CovarianceFull CovarianceMode = "full"
	CovarianceNone CovarianceMode = "none"
)

// Options controls numerical validation. Zero values select documented defaults.
type Options struct {
	// SymmetryTolerance is the absolute tolerance used when checking covariance
	// blocks. The default is 1e-12.
	SymmetryTolerance float64
	// MinVariance rejects smaller independent variances. The default is 1e-20.
	MinVariance float64
	// Solver selects dense or sparse-PCG solution. The default is dense.
	Solver SolverMethod
	// DenseThreshold is the maximum parameter count solved densely by
	// SolverAuto. The default is 300.
	DenseThreshold int
	// Preconditioner selects scalar Jacobi, consecutive block Jacobi, or IC(0)
	// for PCG. Direct batch use defaults to scalar Jacobi.
	Preconditioner PreconditionerMethod
	// PreconditionerBlockSize is required by block Jacobi. ENU networks set it
	// to three so each free station is one block.
	PreconditionerBlockSize int
	// PreconditionerShift is the relative diagonal stabilization used by IC(0).
	// The default is 1e-9.
	PreconditionerShift float64
	// Covariance controls complete covariance materialization. The default is
	// full to preserve the historical Solve contract.
	Covariance CovarianceMode
	// MaxIterations bounds PCG iterations. Zero selects max(100, 10*n).
	MaxIterations int
	// RelativeTolerance is the PCG relative residual tolerance. The default is
	// 1e-10.
	RelativeTolerance float64
	// AbsoluteTolerance is the PCG absolute residual tolerance. The default is
	// 1e-12.
	AbsoluteTolerance float64
}

func (o Options) withDefaults() Options {
	if o.SymmetryTolerance <= 0 {
		o.SymmetryTolerance = 1e-12
	}
	if o.MinVariance <= 0 {
		o.MinVariance = 1e-20
	}
	if o.Solver == "" {
		o.Solver = SolverDense
	}
	if o.DenseThreshold == 0 {
		o.DenseThreshold = 300
	}
	if o.Preconditioner == "" {
		o.Preconditioner = PreconditionerJacobi
	}
	if o.Covariance == "" {
		o.Covariance = CovarianceFull
	}
	if o.Preconditioner == PreconditionerIC0 && o.PreconditionerShift == 0 {
		o.PreconditionerShift = 1e-9
	}
	if o.RelativeTolerance == 0 {
		o.RelativeTolerance = 1e-10
	}
	if o.AbsoluteTolerance == 0 {
		o.AbsoluteTolerance = 1e-12
	}
	return o
}

func normalizeOptions(options *Options) (Options, error) {
	normalized := Options{}
	if options != nil {
		normalized = *options
	}
	if normalized.SymmetryTolerance < 0 || math.IsNaN(normalized.SymmetryTolerance) || math.IsInf(normalized.SymmetryTolerance, 0) {
		return Options{}, invalid("symmetry tolerance must be finite and non-negative")
	}
	if normalized.MinVariance < 0 || math.IsNaN(normalized.MinVariance) || math.IsInf(normalized.MinVariance, 0) {
		return Options{}, invalid("minimum variance must be finite and non-negative")
	}
	if normalized.Solver != "" && normalized.Solver != SolverDense && normalized.Solver != SolverAuto && normalized.Solver != SolverPCG {
		return Options{}, invalid("unsupported solver %q", normalized.Solver)
	}
	if normalized.DenseThreshold < 0 {
		return Options{}, invalid("dense threshold must be non-negative")
	}
	if normalized.Preconditioner != "" && normalized.Preconditioner != PreconditionerJacobi && normalized.Preconditioner != PreconditionerBlockJacobi && normalized.Preconditioner != PreconditionerIC0 {
		return Options{}, invalid("unsupported preconditioner %q", normalized.Preconditioner)
	}
	if normalized.PreconditionerBlockSize < 0 {
		return Options{}, invalid("preconditioner block size must be non-negative")
	}
	if normalized.Preconditioner == PreconditionerBlockJacobi && normalized.PreconditionerBlockSize == 0 {
		return Options{}, invalid("block-Jacobi preconditioner requires a positive block size")
	}
	if normalized.PreconditionerShift < 0 || math.IsNaN(normalized.PreconditionerShift) || math.IsInf(normalized.PreconditionerShift, 0) {
		return Options{}, invalid("preconditioner shift must be finite and non-negative")
	}
	if normalized.Covariance != "" && normalized.Covariance != CovarianceFull && normalized.Covariance != CovarianceNone {
		return Options{}, invalid("unsupported covariance mode %q", normalized.Covariance)
	}
	if normalized.MaxIterations < 0 {
		return Options{}, invalid("maximum iterations must be non-negative")
	}
	if normalized.RelativeTolerance < 0 || math.IsNaN(normalized.RelativeTolerance) || math.IsInf(normalized.RelativeTolerance, 0) {
		return Options{}, invalid("relative tolerance must be finite and non-negative")
	}
	if normalized.AbsoluteTolerance < 0 || math.IsNaN(normalized.AbsoluteTolerance) || math.IsInf(normalized.AbsoluteTolerance, 0) {
		return Options{}, invalid("absolute tolerance must be finite and non-negative")
	}
	return normalized.withDefaults(), nil
}
