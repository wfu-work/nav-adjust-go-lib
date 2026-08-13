package core

import "gonum.org/v1/gonum/mat"

// Residual contains diagnostics for one input observation.
type Residual struct {
	ID           string
	Group        string
	Value        float64
	Variance     float64
	Redundancy   float64
	Standardized float64
}

// ConstraintResidual contains diagnostics for one soft constraint. Exact
// constraints are satisfied by construction and are not included here.
type ConstraintResidual struct {
	ID           string
	Value        float64
	Variance     float64
	Standardized float64
}

// Result contains the parameter correction and adjustment quality information.
//
// FormalCovariance is the unscaled cofactor matrix when CovarianceAvailable is
// true. Covariance equals Sigma0^2*FormalCovariance. Compact solve modes leave
// both matrices nil while retaining raw residuals and global statistics.
type Result struct {
	Delta                        []float64
	FormalCovariance             *mat.SymDense
	Covariance                   *mat.SymDense
	Residuals                    []Residual
	ConstraintResiduals          []ConstraintResidual
	Rank                         int
	DOF                          int
	Sigma0                       float64
	Objective                    float64
	Condition                    float64
	ConditionAvailable           bool
	Method                       string
	SolverPreconditioner         string
	SolverIterations             int
	SolverRelativeResidual       float64
	CovarianceAvailable          bool
	ResidualDiagnosticsAvailable bool
	ObservationCount             int
	SoftConstraintCount          int
	ExactConstraintCount         int
}
