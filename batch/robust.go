package batch

import (
	"context"
	"fmt"
	"math"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
)

// HuberOptions controls generic observation-wise Huber reweighting.
type HuberOptions struct {
	K             float64
	MaxIterations int
	Tolerance     float64
	MinWeight     float64
	Batch         *Options
}

func (options HuberOptions) withDefaults() HuberOptions {
	if options.K <= 0 {
		options.K = 1.5
	}
	if options.MaxIterations <= 0 {
		options.MaxIterations = 10
	}
	if options.Tolerance <= 0 {
		options.Tolerance = 1e-3
	}
	if options.MinWeight <= 0 || options.MinWeight > 1 {
		options.MinWeight = 0.05
	}
	return options
}

func normalizeHuberOptions(options *HuberOptions) (HuberOptions, error) {
	normalized := HuberOptions{}
	if options != nil {
		normalized = *options
	}
	if normalized.K < 0 || math.IsNaN(normalized.K) || math.IsInf(normalized.K, 0) {
		return HuberOptions{}, fmt.Errorf("batch: %w: Huber K must be finite and non-negative", adjust.ErrInvalidProblem)
	}
	if normalized.MaxIterations < 0 {
		return HuberOptions{}, fmt.Errorf("batch: %w: Huber max iterations must be non-negative", adjust.ErrInvalidProblem)
	}
	if normalized.Tolerance < 0 || math.IsNaN(normalized.Tolerance) || math.IsInf(normalized.Tolerance, 0) {
		return HuberOptions{}, fmt.Errorf("batch: %w: Huber tolerance must be finite and non-negative", adjust.ErrInvalidProblem)
	}
	if normalized.MinWeight < 0 || normalized.MinWeight > 1 || math.IsNaN(normalized.MinWeight) || math.IsInf(normalized.MinWeight, 0) {
		return HuberOptions{}, fmt.Errorf("batch: %w: Huber minimum weight must be zero or in (0, 1]", adjust.ErrInvalidProblem)
	}
	return normalized.withDefaults(), nil
}

// HuberResult contains the final adjustment and per-observation Huber weights.
type HuberResult struct {
	Adjustment *adjust.Result
	Weights    []float64
	Iterations int
	Converged  bool
}

// SolveHuber solves a linear problem by inflating observation covariance for
// large normalized residuals. Constraints are not robustly reweighted.
func SolveHuber(problem adjust.Problem, options *HuberOptions) (*HuberResult, error) {
	return SolveHuberContext(context.Background(), problem, options)
}

// SolveHuberContext is SolveHuber with cooperative cancellation across every
// reweighting iteration and batch solve.
func SolveHuberContext(ctx context.Context, problem adjust.Problem, options *HuberOptions) (*HuberResult, error) {
	if ctx == nil {
		return nil, fmt.Errorf("batch: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	if len(problem.Equations) == 0 {
		return nil, fmt.Errorf("batch: %w: Huber problem has no observations", adjust.ErrInvalidProblem)
	}
	opts, err := normalizeHuberOptions(options)
	if err != nil {
		return nil, err
	}
	weights := make([]float64, len(problem.Equations))
	for i := range weights {
		weights[i] = 1
	}
	result := &HuberResult{Weights: weights}
	variances, err := observationVariances(problem)
	if err != nil {
		return nil, err
	}

	for iteration := 1; iteration <= opts.MaxIterations; iteration++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		weighted := reweightObservations(problem, weights)
		adjustment, err := SolveContext(ctx, weighted, opts.Batch)
		if err != nil {
			return nil, fmt.Errorf("batch: Huber iteration %d: %w", iteration, err)
		}
		next := make([]float64, len(weights))
		maxChange := 0.0
		for i, residual := range adjustment.Residuals {
			score := math.Abs(residual.Value) / math.Sqrt(variances[i])
			weight := 1.0
			if score > opts.K {
				weight = math.Max(opts.MinWeight, opts.K/score)
			}
			next[i] = weight
			maxChange = math.Max(maxChange, math.Abs(weight-weights[i]))
		}
		result.Adjustment = adjustment
		result.Iterations = iteration
		result.Weights = next
		weights = next
		if maxChange <= opts.Tolerance {
			result.Converged = true
			if maxChange == 0 {
				return result, nil
			}
			adjustment, err = SolveContext(ctx, reweightObservations(problem, weights), opts.Batch)
			if err != nil {
				return nil, fmt.Errorf("batch: Huber converged solve: %w", err)
			}
			result.Adjustment = adjustment
			return result, nil
		}
	}
	adjustment, err := SolveContext(ctx, reweightObservations(problem, weights), opts.Batch)
	if err != nil {
		return nil, fmt.Errorf("batch: Huber final solve: %w", err)
	}
	result.Adjustment = adjustment
	return result, nil
}

func observationVariances(problem adjust.Problem) ([]float64, error) {
	variances := make([]float64, len(problem.Equations))
	for i, equation := range problem.Equations {
		variances[i] = equation.Variance
	}
	for _, block := range problem.CovarianceBlocks {
		n := len(block.RowIndexes)
		if len(block.Covariance) != n*n {
			return nil, fmt.Errorf("batch: %w: covariance block %q has invalid dimensions", adjust.ErrInvalidProblem, block.ID)
		}
		for local, row := range block.RowIndexes {
			if row < 0 || row >= len(variances) {
				return nil, fmt.Errorf("batch: %w: covariance block %q has an invalid row index", adjust.ErrInvalidProblem, block.ID)
			}
			variances[row] = block.Covariance[local*n+local]
		}
	}
	for i, variance := range variances {
		if variance <= 0 || math.IsNaN(variance) || math.IsInf(variance, 0) {
			return nil, fmt.Errorf("batch: %w: observation %d has invalid variance", adjust.ErrInvalidProblem, i)
		}
	}
	return variances, nil
}

func reweightObservations(problem adjust.Problem, weights []float64) adjust.Problem {
	weighted := adjust.CloneProblem(problem)
	blocked := make([]bool, len(weighted.Equations))
	for blockIndex := range weighted.CovarianceBlocks {
		block := &weighted.CovarianceBlocks[blockIndex]
		n := len(block.RowIndexes)
		for r, rowR := range block.RowIndexes {
			blocked[rowR] = true
			for c, rowC := range block.RowIndexes {
				block.Covariance[r*n+c] /= math.Sqrt(weights[rowR] * weights[rowC])
			}
		}
	}
	for i := range weighted.Equations {
		if !blocked[i] {
			weighted.Equations[i].Variance /= weights[i]
		}
	}
	return weighted
}
