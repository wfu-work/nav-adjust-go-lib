package model

// RobustMethod identifies a robust weighting function.
type RobustMethod string

// SolverMethod selects the numerical normal-equation solver.
type SolverMethod string

// PreconditionerMethod selects the PCG preconditioner.
type PreconditionerMethod string

// DatumMode controls how translational datum defects are handled.
type DatumMode string

const (
	// SolverDense uses dense direct factorization and is the compatibility
	// default for small and medium networks.
	SolverDense SolverMethod = "dense"
	// SolverAuto keeps full-covariance and small systems dense, then switches
	// larger reduced-covariance systems to sparse PCG.
	SolverAuto SolverMethod = "auto"
	// SolverPCG uses a sparse normal matrix and preconditioned conjugate
	// gradients. Exact constraints are handled by null-space projection.
	SolverPCG SolverMethod = "pcg"
)

const (
	// PreconditionerJacobi uses one scalar normal-matrix diagonal entry per
	// parameter.
	PreconditionerJacobi PreconditionerMethod = "jacobi"
	// PreconditionerBlockJacobi uses each station's 3-by-3 ENU diagonal block.
	PreconditionerBlockJacobi PreconditionerMethod = "block-jacobi"
	// PreconditionerIC0 uses an incomplete Cholesky factor with the normal
	// matrix's existing sparse pattern.
	PreconditionerIC0 PreconditionerMethod = "ic0"
)

const (
	// DatumExternal requires every connected component to contain a fixed
	// station or a stochastic position prior. It is the compatibility default.
	DatumExternal DatumMode = "external"
	// DatumFreeCentroid permits components made only from relative baselines and
	// defines their internal datum by constraining the component centroid to
	// zero in East, North, and Up.
	DatumFreeCentroid DatumMode = "free-centroid"
)

// CovarianceMode selects how much covariance information is returned.
type CovarianceMode string

const (
	// CovarianceFull returns the complete parameter covariance, station blocks,
	// and observation residual diagnostics. It is the compatibility default.
	CovarianceFull CovarianceMode = "full"
	// CovarianceStationBlocks omits the complete n-by-n matrices but computes
	// station 3-by-3 blocks and observation residual diagnostics on demand.
	CovarianceStationBlocks CovarianceMode = "station-blocks"
	// CovarianceNone returns coordinates, raw residuals, objective, variance
	// factor, and global test without covariance-derived diagnostics.
	CovarianceNone CovarianceMode = "none"
)

const (
	// RobustHuber selects block-wise Huber weighting of complete ENU baselines.
	RobustHuber RobustMethod = "huber"
)

// RobustOptions controls optional baseline-wise robust reweighting. Zero
// values select documented defaults.
type RobustOptions struct {
	Method        RobustMethod `json:"method,omitempty"`
	Threshold     float64      `json:"threshold,omitempty"`
	MaxIterations int          `json:"max_iterations,omitempty"`
	Tolerance     float64      `json:"tolerance,omitempty"`
	MinWeight     float64      `json:"min_weight,omitempty"`
}

// SolverOptions controls direct, automatic, and sparse iterative solution.
// Zero values select the dense solver and documented PCG tolerances.
type SolverOptions struct {
	Method              SolverMethod         `json:"method,omitempty"`
	DenseThreshold      int                  `json:"dense_threshold,omitempty"`
	Preconditioner      PreconditionerMethod `json:"preconditioner,omitempty"`
	PreconditionerShift float64              `json:"preconditioner_shift,omitempty"`
	MaxIterations       int                  `json:"max_iterations,omitempty"`
	RelativeTolerance   float64              `json:"relative_tolerance,omitempty"`
	AbsoluteTolerance   float64              `json:"absolute_tolerance,omitempty"`
}

// VarianceComponentOptions controls optional baseline-group covariance-scale
// estimation. Each distinct Baseline.Group is treated as one component; an
// empty group name is a valid default component. Zero values select documented
// defaults.
type VarianceComponentOptions struct {
	MaxIterations     int     `json:"max_iterations,omitempty"`
	Tolerance         float64 `json:"tolerance,omitempty"`
	MinScale          float64 `json:"min_scale,omitempty"`
	MaxScale          float64 `json:"max_scale,omitempty"`
	MinimumRedundancy float64 `json:"minimum_redundancy,omitempty"`
}

// ENUNetworkOptions controls datum definition, numerical solution, covariance
// output, and optional robust or variance-component estimation.
type ENUNetworkOptions struct {
	Robust             *RobustOptions            `json:"robust,omitempty"`
	VarianceComponents *VarianceComponentOptions `json:"variance_components,omitempty"`
	Solver             SolverOptions             `json:"solver,omitempty"`
	Covariance         CovarianceMode            `json:"covariance,omitempty"`
	Datum              DatumMode                 `json:"datum,omitempty"`
	Confidence         float64                   `json:"confidence,omitempty"`
	SymmetryTolerance  float64                   `json:"symmetry_tolerance,omitempty"`
	MinimumVariance    float64                   `json:"minimum_variance,omitempty"`
	ConditionWarnLimit float64                   `json:"condition_warn_limit,omitempty"`
}

// Options is the short name used with network.Solve and network.Validate.
type Options = ENUNetworkOptions
