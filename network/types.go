// Package network validates and adjusts networks of relative ENU vectors.
package network

import "github.com/wfu-work/nav-adjust-go-lib/model"

// Public network types are aliases of model contracts. Keeping the contracts
// in model lets applications construct or serialize data without importing the
// numerical solver package.
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
	Baseline                 = model.Baseline
	Prior                    = model.Prior
	Problem                  = model.Problem
	Options                  = model.Options
	Result                   = model.Result
)

const (
	// RobustHuber selects block-wise Huber weighting of complete ENU baselines.
	RobustHuber               = model.RobustHuber
	SolverDense               = model.SolverDense
	SolverPCG                 = model.SolverPCG
	PreconditionerJacobi      = model.PreconditionerJacobi
	PreconditionerBlockJacobi = model.PreconditionerBlockJacobi
	CovarianceFull            = model.CovarianceFull
	CovarianceStationBlocks   = model.CovarianceStationBlocks
	CovarianceNone            = model.CovarianceNone
	DatumExternal             = model.DatumExternal
	DatumFreeCentroid         = model.DatumFreeCentroid
)

// DiagonalMatrix3 creates a diagonal 3-by-3 covariance from variances.
func DiagonalMatrix3(eastVariance, northVariance, upVariance float64) Matrix3 {
	return model.DiagonalMatrix3(eastVariance, northVariance, upVariance)
}

// Matrix3FromStdDev creates a diagonal covariance from standard deviations.
func Matrix3FromStdDev(east, north, up float64) Matrix3 {
	return model.Matrix3FromStdDev(east, north, up)
}
