package network

import (
	"context"
	"fmt"
	"math"
	"sort"

	"github.com/wfu-work/nav-adjust-go-lib/batch"
	"github.com/wfu-work/nav-adjust-go-lib/core"
	"github.com/wfu-work/nav-adjust-go-lib/internal/statistics"
	"gonum.org/v1/gonum/mat"
)

type networkSolve struct {
	adjustment         *core.Result
	detailed           *batch.Solution
	weights            []float64
	iterations         int
	converged          bool
	problem            core.Problem
	varianceComponents []VarianceComponentResult
	varianceIterations int
	varianceConverged  bool
}

// SolveENUNetwork adjusts relative From -> To ENU vectors into station
// positions in one common ENU frame.
func SolveENUNetwork(input ENUNetworkProblem, options *ENUNetworkOptions) (*ENUNetworkResult, error) {
	return SolveENUNetworkContext(context.Background(), input, options)
}

// SolveENUNetworkContext adjusts an ENU network and supports cancellation of
// all iterative solves and on-demand covariance queries.
func SolveENUNetworkContext(ctx context.Context, input ENUNetworkProblem, options *ENUNetworkOptions) (*ENUNetworkResult, error) {
	if ctx == nil {
		return nil, fmt.Errorf("adjust: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	opts, err := normalizeNetworkOptions(options)
	if err != nil {
		return nil, err
	}
	compiled, err := compileENUNetwork(input, opts)
	if err != nil {
		return nil, err
	}
	solution, err := solveCompiledNetwork(ctx, compiled, len(input.Baselines))
	if err != nil {
		return nil, err
	}
	return buildNetworkResult(ctx, input, compiled, solution)
}

// Solve adjusts relative From -> To ENU vectors into station positions.
func Solve(input Problem, options *Options) (*Result, error) {
	return SolveENUNetwork(input, options)
}

// SolveContext is Solve with cooperative cancellation.
func SolveContext(ctx context.Context, input Problem, options *Options) (*Result, error) {
	return SolveENUNetworkContext(ctx, input, options)
}

func solveCompiledNetwork(ctx context.Context, compiled compiledNetwork, baselineCount int) (networkSolve, error) {
	if compiled.publicOptions.variance != nil && baselineCount > 0 {
		return solveVarianceComponentNetwork(ctx, compiled, baselineCount)
	}
	return solveNetworkProblem(ctx, compiled, compiled.problem, baselineCount, compiled.publicOptions.batch)
}

func solveNetworkProblem(ctx context.Context, compiled compiledNetwork, problem core.Problem, baselineCount int, finalOptions batch.Options) (networkSolve, error) {
	if compiled.publicOptions.robust == nil || baselineCount == 0 {
		detailed, err := batch.SolveDetailedContext(ctx, problem, &finalOptions)
		if err != nil {
			return networkSolve{}, fmt.Errorf("adjust: solve ENU network: %w", err)
		}
		weights := make([]float64, baselineCount)
		for i := range weights {
			weights[i] = 1
		}
		return networkSolve{adjustment: detailed.Result, detailed: detailed, weights: weights, iterations: 1, converged: true, problem: problem}, nil
	}
	return solveRobustNetwork(ctx, compiled, problem, baselineCount, finalOptions)
}

func solveRobustNetwork(ctx context.Context, compiled compiledNetwork, baseProblem core.Problem, baselineCount int, finalOptions batch.Options) (networkSolve, error) {
	options := compiled.publicOptions.robust
	weights := make([]float64, baselineCount)
	for i := range weights {
		weights[i] = 1
	}
	result := networkSolve{weights: weights}
	iterationOptions := finalOptions
	iterationOptions.Covariance = batch.CovarianceNone
	for iteration := 1; iteration <= options.maxIterations; iteration++ {
		if err := ctx.Err(); err != nil {
			return networkSolve{}, err
		}
		weighted := applyBaselineWeights(baseProblem, weights)
		adjustment, err := batch.SolveContext(ctx, weighted, &iterationOptions)
		if err != nil {
			return networkSolve{}, fmt.Errorf("adjust: robust ENU iteration %d: %w", iteration, err)
		}
		next := make([]float64, baselineCount)
		maxChange := 0.0
		for baselineIndex := range baselineCount {
			if err := ctx.Err(); err != nil {
				return networkSolve{}, err
			}
			residual := []float64{
				adjustment.Residuals[baselineIndex*3].Value,
				adjustment.Residuals[baselineIndex*3+1].Value,
				adjustment.Residuals[baselineIndex*3+2].Value,
			}
			score, err := mahalanobis3(residual, baseProblem.CovarianceBlocks[baselineIndex].Covariance)
			if err != nil {
				return networkSolve{}, err
			}
			weight := 1.0
			if score > options.threshold {
				weight = math.Max(options.minWeight, options.threshold/score)
			}
			next[baselineIndex] = weight
			maxChange = math.Max(maxChange, math.Abs(weight-weights[baselineIndex]))
		}
		weights = next
		result = networkSolve{adjustment: adjustment, weights: weights, iterations: iteration, problem: weighted}
		if maxChange <= options.tolerance {
			result.converged = true
			result.problem = applyBaselineWeights(baseProblem, weights)
			result.detailed, err = batch.SolveDetailedContext(ctx, result.problem, &finalOptions)
			if err != nil {
				return networkSolve{}, fmt.Errorf("adjust: robust ENU converged solve: %w", err)
			}
			result.adjustment = result.detailed.Result
			return result, nil
		}
	}
	result.problem = applyBaselineWeights(baseProblem, weights)
	detailed, err := batch.SolveDetailedContext(ctx, result.problem, &finalOptions)
	if err != nil {
		return networkSolve{}, fmt.Errorf("adjust: robust ENU final solve: %w", err)
	}
	result.detailed = detailed
	result.adjustment = detailed.Result
	return result, nil
}

func solveVarianceComponentNetwork(ctx context.Context, compiled compiledNetwork, baselineCount int) (networkSolve, error) {
	options := compiled.publicOptions.variance
	groups := baselineGroups(compiled.problem, baselineCount)
	scales := make(map[string]float64, len(groups))
	for _, group := range groups {
		scales[group] = 1
	}
	varianceOptions := statistics.VarianceOptions{
		Tolerance: options.tolerance, MinScale: options.minScale, MaxScale: options.maxScale,
		MinimumRedundancy: options.minimumRedundancy,
	}
	iterations := 0
	converged := false
	for iteration := 1; iteration <= options.maxIterations; iteration++ {
		if err := ctx.Err(); err != nil {
			return networkSolve{}, err
		}
		iterations = iteration
		scaled := applyBaselineGroupScales(compiled.problem, scales, baselineCount)
		iterationOptions := compiled.publicOptions.batch
		iterationOptions.Covariance = batch.CovarianceNone
		solution, err := solveNetworkProblem(ctx, compiled, scaled, baselineCount, iterationOptions)
		if err != nil {
			return networkSolve{}, fmt.Errorf("adjust: variance-component iteration %d: %w", iteration, err)
		}
		components, err := varianceComponentStatistics(ctx, compiled, solution, groups, scales, baselineCount)
		if err != nil {
			return networkSolve{}, fmt.Errorf("adjust: variance-component iteration %d: %w", iteration, err)
		}
		update, err := statistics.UpdateGroupScales(components, varianceOptions)
		if err != nil {
			return networkSolve{}, fmt.Errorf("adjust: variance-component iteration %d: %w", iteration, err)
		}
		for _, component := range update.Components {
			scales[component.ID] = component.Scale
		}
		if update.Converged {
			converged = true
			break
		}
	}

	finalProblem := applyBaselineGroupScales(compiled.problem, scales, baselineCount)
	result, err := solveNetworkProblem(ctx, compiled, finalProblem, baselineCount, compiled.publicOptions.batch)
	if err != nil {
		return networkSolve{}, fmt.Errorf("adjust: final variance-component solve: %w", err)
	}
	components, err := varianceComponentStatistics(ctx, compiled, result, groups, scales, baselineCount)
	if err != nil {
		return networkSolve{}, fmt.Errorf("adjust: final variance-component statistics: %w", err)
	}
	result.varianceComponents = publicVarianceComponents(components)
	result.varianceIterations = iterations
	result.varianceConverged = converged
	return result, nil
}

func baselineGroups(problem core.Problem, baselineCount int) []string {
	set := make(map[string]struct{})
	for i := 0; i < baselineCount; i++ {
		block := problem.CovarianceBlocks[i]
		set[problem.Equations[block.RowIndexes[0]].Group] = struct{}{}
	}
	groups := make([]string, 0, len(set))
	for group := range set {
		groups = append(groups, group)
	}
	sort.Strings(groups)
	return groups
}

func applyBaselineGroupScales(problem core.Problem, scales map[string]float64, baselineCount int) core.Problem {
	scaled := core.CloneProblem(problem)
	for i := 0; i < baselineCount; i++ {
		block := &scaled.CovarianceBlocks[i]
		group := scaled.Equations[block.RowIndexes[0]].Group
		for j := range block.Covariance {
			block.Covariance[j] *= scales[group]
		}
	}
	return scaled
}

func varianceComponentStatistics(
	ctx context.Context,
	compiled compiledNetwork,
	solution networkSolve,
	groups []string,
	scales map[string]float64,
	baselineCount int,
) ([]statistics.VarianceComponent, error) {
	formal, err := queryNetworkCovariance(ctx, solution.detailed, solution.problem, compiled.stationIndex)
	if err != nil {
		return nil, err
	}
	byGroup := make(map[string]*statistics.VarianceComponent, len(groups))
	for _, group := range groups {
		byGroup[group] = &statistics.VarianceComponent{ID: group, Scale: scales[group]}
	}
	for i := 0; i < baselineCount; i++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		block := solution.problem.CovarianceBlocks[i]
		group := solution.problem.Equations[block.RowIndexes[0]].Group
		component := byGroup[group]
		component.BaselineCount++
		component.ObservationCount += len(block.RowIndexes)
		residual := []float64{
			solution.adjustment.Residuals[i*3].Value,
			solution.adjustment.Residuals[i*3+1].Value,
			solution.adjustment.Residuals[i*3+2].Value,
		}
		formalResidual := observationResidualCovariance(solution.problem, formal, i)
		objective, redundancy, err := blockVarianceStatistics(residual, block.Covariance, formalResidual)
		if err != nil {
			return nil, err
		}
		component.Objective += objective
		component.Redundancy += redundancy
	}
	components := make([]statistics.VarianceComponent, len(groups))
	for i, group := range groups {
		components[i] = *byGroup[group]
	}
	return components, nil
}

func blockVarianceStatistics(residual, covariance []float64, residualCovariance Matrix3) (float64, float64, error) {
	observation := mat.NewSymDense(3, nil)
	formalResidual := mat.NewSymDense(3, nil)
	for row := range 3 {
		for column := 0; column <= row; column++ {
			observation.SetSym(row, column, covariance[row*3+column])
			formalResidual.SetSym(row, column, residualCovariance.At(row, column))
		}
	}
	var factor mat.Cholesky
	if !factor.Factorize(observation) {
		return 0, 0, fmt.Errorf("%w: variance-component covariance", ErrNotPositiveDefinite)
	}
	residualVector := mat.NewVecDense(3, residual)
	weightedResidual := mat.NewVecDense(3, nil)
	if err := factor.SolveVecTo(weightedResidual, residualVector); err != nil {
		return 0, 0, fmt.Errorf("adjust: variance-component objective: %w", err)
	}
	var weightedCovariance mat.Dense
	if err := factor.SolveTo(&weightedCovariance, formalResidual); err != nil {
		return 0, 0, fmt.Errorf("adjust: variance-component redundancy: %w", err)
	}
	redundancy := 0.0
	for i := range 3 {
		redundancy += weightedCovariance.At(i, i)
	}
	return math.Max(0, mat.Dot(residualVector, weightedResidual)), math.Max(0, redundancy), nil
}

func publicVarianceComponents(components []statistics.VarianceComponent) []VarianceComponentResult {
	result := make([]VarianceComponentResult, len(components))
	for i, component := range components {
		result[i] = VarianceComponentResult{
			Group: component.ID, Scale: component.Scale, StdDevScale: math.Sqrt(component.Scale),
			BaselineCount: component.BaselineCount, ObservationCount: component.ObservationCount,
			Objective: component.Objective, Redundancy: component.Redundancy,
		}
	}
	return result
}

func applyBaselineWeights(problem core.Problem, weights []float64) core.Problem {
	weighted := core.CloneProblem(problem)
	for i := range weights {
		for j := range weighted.CovarianceBlocks[i].Covariance {
			weighted.CovarianceBlocks[i].Covariance[j] /= weights[i]
		}
	}
	return weighted
}

func mahalanobis3(residual, covariance []float64) (float64, error) {
	matrix := mat.NewSymDense(3, nil)
	for row := range 3 {
		for column := 0; column <= row; column++ {
			matrix.SetSym(row, column, covariance[row*3+column])
		}
	}
	var chol mat.Cholesky
	if !chol.Factorize(matrix) {
		return 0, fmt.Errorf("%w: baseline covariance", ErrNotPositiveDefinite)
	}
	vector := mat.NewVecDense(3, residual)
	weighted := mat.NewVecDense(3, nil)
	if err := chol.SolveVecTo(weighted, vector); err != nil {
		return 0, fmt.Errorf("adjust: robust residual score: %w", err)
	}
	return math.Sqrt(math.Max(0, mat.Dot(vector, weighted))), nil
}
