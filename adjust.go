// Package adjust is the stable facade of nav-adjust-go-lib. New applications
// can use this root package directly; model and network provide separated
// contracts and solver layers for larger projects.
package adjust

import (
	"context"

	"github.com/wfu-work/nav-adjust-go-lib/model"
	"github.com/wfu-work/nav-adjust-go-lib/network"
)

// Root types alias the standalone model contracts, so callers can move between
// the facade and subpackages without conversion or allocation.
type (
	ENU                      = model.ENU
	Matrix                   = model.Matrix
	Matrix3                  = model.Matrix3
	Station                  = model.Station
	ENUBaseline              = model.ENUBaseline
	PositionPrior            = model.PositionPrior
	ENUNetworkProblem        = model.ENUNetworkProblem
	RobustMethod             = model.RobustMethod
	RobustOptions            = model.RobustOptions
	SolverMethod             = model.SolverMethod
	SolverOptions            = model.SolverOptions
	PreconditionerMethod     = model.PreconditionerMethod
	DatumMode                = model.DatumMode
	VarianceComponentOptions = model.VarianceComponentOptions
	CovarianceMode           = model.CovarianceMode
	ENUNetworkOptions        = model.ENUNetworkOptions
	AdjustedStation          = model.AdjustedStation
	BaselineResult           = model.BaselineResult
	PositionPriorResult      = model.PositionPriorResult
	GlobalTestResult         = model.GlobalTestResult
	VarianceComponentResult  = model.VarianceComponentResult
	NetworkDiagnostics       = model.NetworkDiagnostics
	ENUNetworkResult         = model.ENUNetworkResult
	ValidationError          = network.ValidationError
)

const (
	// RobustHuber selects block-wise Huber weighting of complete ENU baselines.
	RobustHuber               = model.RobustHuber
	SolverDense               = model.SolverDense
	SolverAuto                = model.SolverAuto
	SolverPCG                 = model.SolverPCG
	PreconditionerJacobi      = model.PreconditionerJacobi
	PreconditionerBlockJacobi = model.PreconditionerBlockJacobi
	PreconditionerIC0         = model.PreconditionerIC0
	CovarianceFull            = model.CovarianceFull
	CovarianceStationBlocks   = model.CovarianceStationBlocks
	CovarianceNone            = model.CovarianceNone
	DatumExternal             = model.DatumExternal
	DatumFreeCentroid         = model.DatumFreeCentroid
)

var (
	ErrInvalidProblem         = network.ErrInvalidProblem
	ErrRankDeficient          = network.ErrRankDeficient
	ErrNotPositiveDefinite    = network.ErrNotPositiveDefinite
	ErrInvalidCovariance      = network.ErrInvalidCovariance
	ErrDuplicateStation       = network.ErrDuplicateStation
	ErrDuplicateBaseline      = network.ErrDuplicateBaseline
	ErrDuplicatePrior         = network.ErrDuplicatePrior
	ErrUnknownStation         = network.ErrUnknownStation
	ErrDisconnectedNetwork    = network.ErrDisconnectedNetwork
	ErrUnsupportedMethod      = network.ErrUnsupportedMethod
	ErrNotConverged           = network.ErrNotConverged
	ErrInsufficientRedundancy = network.ErrInsufficientRedundancy
)

// DiagonalMatrix3 creates a diagonal 3-by-3 covariance from variances.
func DiagonalMatrix3(eastVariance, northVariance, upVariance float64) Matrix3 {
	return model.DiagonalMatrix3(eastVariance, northVariance, upVariance)
}

// Matrix3FromStdDev creates a diagonal covariance from standard deviations.
func Matrix3FromStdDev(east, north, up float64) Matrix3 {
	return model.Matrix3FromStdDev(east, north, up)
}

// ValidateENUNetwork validates input without running the adjustment.
func ValidateENUNetwork(problem ENUNetworkProblem, options *ENUNetworkOptions) error {
	return network.ValidateENUNetwork(problem, options)
}

// SolveENUNetwork adjusts relative ENU vectors into station coordinates.
func SolveENUNetwork(problem ENUNetworkProblem, options *ENUNetworkOptions) (*ENUNetworkResult, error) {
	return network.SolveENUNetwork(problem, options)
}

// SolveENUNetworkContext adjusts relative ENU vectors and permits callers to
// cancel sparse iteration, robust reweighting, covariance queries, and
// variance-component estimation.
func SolveENUNetworkContext(ctx context.Context, problem ENUNetworkProblem, options *ENUNetworkOptions) (*ENUNetworkResult, error) {
	return network.SolveENUNetworkContext(ctx, problem, options)
}
