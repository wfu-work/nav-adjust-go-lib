package network

import (
	"context"
	"fmt"
	"math"

	"github.com/wfu-work/nav-adjust-go-lib/quality"
)

func buildNetworkResult(ctx context.Context, input ENUNetworkProblem, compiled compiledNetwork, solution networkSolve) (*ENUNetworkResult, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	adjustment := solution.adjustment
	positions := make(map[string]ENU, len(input.Stations))
	covarianceMode := compiled.publicOptions.covariance
	var formalCovariance covarianceReader
	switch covarianceMode {
	case CovarianceFull:
		formalCovariance = adjustment.FormalCovariance
	case CovarianceStationBlocks:
		queried, err := queryNetworkCovariance(ctx, solution.detailed, solution.problem, compiled.stationIndex)
		if err != nil {
			return nil, err
		}
		formalCovariance = queried
	}
	localCovarianceAvailable := formalCovariance != nil
	result := &ENUNetworkResult{
		Name:          input.Name,
		ParameterKeys: append([]string(nil), compiled.parameterKeys...),
		Diagnostics: NetworkDiagnostics{
			StationCount: len(input.Stations), BaselineCount: len(input.Baselines), PriorCount: len(input.Priors),
			ObservationCount: (len(input.Baselines) + len(input.Priors)) * 3, ParameterCount: len(compiled.parameterKeys),
			Rank: min(adjustment.Rank, len(compiled.parameterKeys)), DegreesOfFreedom: adjustment.DOF,
			Objective: adjustment.Objective, Sigma0: adjustment.Sigma0,
			ConditionNumber: adjustment.Condition, Solver: adjustment.Method,
			SolverPreconditioner:         adjustment.SolverPreconditioner,
			ConditionNumberAvailable:     adjustment.ConditionAvailable,
			SolverIterations:             adjustment.SolverIterations,
			SolverRelativeResidual:       adjustment.SolverRelativeResidual,
			CovarianceMode:               covarianceMode,
			FullCovarianceAvailable:      covarianceMode == CovarianceFull,
			StationCovarianceAvailable:   localCovarianceAvailable,
			ResidualDiagnosticsAvailable: localCovarianceAvailable,
			DatumMode:                    compiled.publicOptions.datum,
			FreeDatumComponentCount:      len(compiled.freeDatumComponents),
			InternalDatumConstraintCount: len(compiled.freeDatumComponents) * 3,
			Iterations:                   solution.iterations, Converged: solution.converged,
			VarianceComponentIterations: solution.varianceIterations,
			VarianceComponentsAvailable: compiled.publicOptions.variance != nil,
			VarianceComponentsConverged: solution.varianceConverged,
		},
		VarianceComponents: append([]VarianceComponentResult(nil), solution.varianceComponents...),
	}
	varianceScales := make(map[string]float64, len(solution.varianceComponents))
	for _, component := range solution.varianceComponents {
		varianceScales[component.Group] = component.Scale
	}
	if covarianceMode == CovarianceFull && len(compiled.parameterKeys) > 0 {
		result.FormalCovariance = publicMatrix(adjustment.FormalCovariance, len(compiled.parameterKeys))
		result.Covariance = publicMatrix(adjustment.Covariance, len(compiled.parameterKeys))
	} else {
		result.FormalCovariance = Matrix{Rows: 0, Cols: 0, Data: []float64{}}
		result.Covariance = Matrix{Rows: 0, Cols: 0, Data: []float64{}}
	}
	for _, station := range input.Stations {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		position := stationPosition(station)
		stationResult := AdjustedStation{
			ID: station.ID, Name: station.Name, Fixed: station.Fixed,
			Position: position, Metadata: cloneMetadata(station.Metadata),
		}
		if station.Fixed {
			result.Diagnostics.FixedStationCount++
		} else {
			result.Diagnostics.FreeStationCount++
			indexes := compiled.stationIndex[station.ID]
			position = ENU{
				East: adjustment.Delta[indexes[0]], North: adjustment.Delta[indexes[1]], Up: adjustment.Delta[indexes[2]],
			}
			stationResult.Position = position
			if localCovarianceAvailable {
				stationResult.FormalCovariance = covarianceBlock3(formalCovariance, indexes)
				stationResult.Covariance = scaleMatrix3(stationResult.FormalCovariance, adjustment.Sigma0*adjustment.Sigma0)
				stationResult.StdDev = ENU{
					East: safeSqrt(stationResult.Covariance.At(0, 0)), North: safeSqrt(stationResult.Covariance.At(1, 1)),
					Up: safeSqrt(stationResult.Covariance.At(2, 2)),
				}
			}
		}
		positions[station.ID] = position
		result.Stations = append(result.Stations, stationResult)
	}

	for i, baseline := range input.Baselines {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		from := positions[baseline.From]
		to := positions[baseline.To]
		adjusted := subtractENU(to, from)
		residual := subtractENU(baseline.Vector, adjusted)
		formalResidualCovariance := Matrix3{}
		residualCovariance := Matrix3{}
		standardized := ENU{}
		redundancy := ENU{}
		if localCovarianceAvailable {
			formalResidualCovariance = baselineResidualCovariance(solution.problem, formalCovariance, i)
			residualCovariance = scaleMatrix3(formalResidualCovariance, adjustment.Sigma0*adjustment.Sigma0)
			var err error
			standardized, redundancy, err = observationDiagnostics(
				formalResidualCovariance, solution.problem.CovarianceBlocks[i].Covariance, residual, adjustment.Sigma0,
			)
			if err != nil {
				return nil, err
			}
		}
		weight := solution.weights[i]
		varianceScale := 1.0
		if scale, exists := varianceScales[baseline.Group]; exists {
			varianceScale = scale
		}
		result.Baselines = append(result.Baselines, BaselineResult{
			ID: baseline.ID, From: baseline.From, To: baseline.To, Group: baseline.Group,
			Observed: baseline.Vector, Adjusted: adjusted, Residual: residual,
			ResidualStdDev: ENU{
				East:  safeSqrt(residualCovariance.At(0, 0)),
				North: safeSqrt(residualCovariance.At(1, 1)),
				Up:    safeSqrt(residualCovariance.At(2, 2)),
			},
			Standardized: standardized, Redundancy: redundancy,
			FormalResidualCovariance: formalResidualCovariance,
			ResidualCovariance:       residualCovariance,
			Weight:                   weight, Downweighted: weight < 1,
			VarianceScale: varianceScale,
			Metadata:      cloneMetadata(baseline.Metadata),
		})
	}

	for i, prior := range input.Priors {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		adjusted := positions[prior.StationID]
		residual := subtractENU(prior.Position, adjusted)
		blockIndex := len(input.Baselines) + i
		formalResidualCovariance := Matrix3{}
		residualCovariance := Matrix3{}
		standardized := ENU{}
		redundancy := ENU{}
		if localCovarianceAvailable {
			formalResidualCovariance = observationResidualCovariance(solution.problem, formalCovariance, blockIndex)
			residualCovariance = scaleMatrix3(formalResidualCovariance, adjustment.Sigma0*adjustment.Sigma0)
			var err error
			standardized, redundancy, err = observationDiagnostics(
				formalResidualCovariance, solution.problem.CovarianceBlocks[blockIndex].Covariance, residual, adjustment.Sigma0,
			)
			if err != nil {
				return nil, err
			}
		}
		result.Priors = append(result.Priors, PositionPriorResult{
			ID: prior.ID, StationID: prior.StationID,
			Observed: prior.Position, Adjusted: adjusted, Residual: residual,
			ResidualStdDev: ENU{
				East:  safeSqrt(residualCovariance.At(0, 0)),
				North: safeSqrt(residualCovariance.At(1, 1)),
				Up:    safeSqrt(residualCovariance.At(2, 2)),
			},
			Standardized: standardized, Redundancy: redundancy,
			FormalResidualCovariance: formalResidualCovariance,
			ResidualCovariance:       residualCovariance,
			Metadata:                 cloneMetadata(prior.Metadata),
		})
	}
	reweighted := hasDownweightedBaseline(solution.weights)
	varianceEstimated := compiled.publicOptions.variance != nil
	if adjustment.DOF > 0 && !reweighted && !varianceEstimated {
		global, err := quality.ChiSquare(adjustment, compiled.publicOptions.confidence)
		if err != nil {
			return nil, err
		}
		result.Diagnostics.GlobalTest = &GlobalTestResult{
			Statistic: global.Statistic, DOF: global.DOF, Confidence: global.Confidence,
			Lower: global.Lower, Upper: global.Upper, Passed: global.Passed,
		}
	}
	if adjustment.DOF > 0 && reweighted {
		result.Warnings = append(result.Warnings, "global chi-square test is omitted after data-dependent robust reweighting")
	}
	if adjustment.DOF > 0 && varianceEstimated {
		result.Warnings = append(result.Warnings, "global chi-square test is omitted after data-dependent variance-component estimation")
	}
	if adjustment.ConditionAvailable && adjustment.Condition > compiled.publicOptions.conditionLimit {
		result.Warnings = append(result.Warnings, fmt.Sprintf("condition number %.6g exceeds warning limit %.6g", adjustment.Condition, compiled.publicOptions.conditionLimit))
	}
	if !solution.converged {
		result.Warnings = append(result.Warnings, "robust adjustment reached the iteration limit")
	}
	if varianceEstimated && !solution.varianceConverged {
		result.Warnings = append(result.Warnings, "variance-component estimation reached the iteration limit")
	}
	if len(compiled.freeDatumComponents) > 0 {
		result.Warnings = append(result.Warnings, fmt.Sprintf(
			"free-network centroid datum applied to %d unanchored component(s); coordinates are internal-datum coordinates",
			len(compiled.freeDatumComponents),
		))
	}
	return result, nil
}

func hasDownweightedBaseline(weights []float64) bool {
	for _, weight := range weights {
		if weight < 1 {
			return true
		}
	}
	return false
}

func subtractENU(left, right ENU) ENU {
	return ENU{East: left.East - right.East, North: left.North - right.North, Up: left.Up - right.Up}
}

func safeSqrt(value float64) float64 {
	if value < 0 && value > -1e-12 {
		value = 0
	}
	return math.Sqrt(math.Max(0, value))
}

func cloneMetadata(src map[string]string) map[string]string {
	if src == nil {
		return nil
	}
	dst := make(map[string]string, len(src))
	for key, value := range src {
		dst[key] = value
	}
	return dst
}
