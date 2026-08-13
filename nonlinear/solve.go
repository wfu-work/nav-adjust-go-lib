// Package nonlinear implements iterative linearization on top of the batch
// least-squares solver.
package nonlinear

import (
	"context"
	"fmt"
	"math"
	"reflect"

	"github.com/wfu-work/nav-adjust-go-lib/batch"
	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
)

// Model builds correction equations around state. Each returned equation must
// follow w = H*dx + e, and its parameter count must equal len(state).
type Model interface {
	Linearize(state []float64) (adjust.Problem, error)
}

// ModelFunc adapts a function to Model.
type ModelFunc func(state []float64) (adjust.Problem, error)

func (f ModelFunc) Linearize(state []float64) (adjust.Problem, error) { return f(state) }

// Options controls Gauss-Newton iteration.
type Options struct {
	MaxIterations      int
	StepTolerance      float64
	ObjectiveTolerance float64
	Batch              *batch.Options
}

func (options Options) withDefaults() Options {
	if options.MaxIterations <= 0 {
		options.MaxIterations = 10
	}
	if options.StepTolerance <= 0 {
		options.StepTolerance = 1e-6
	}
	if options.ObjectiveTolerance <= 0 {
		options.ObjectiveTolerance = 1e-10
	}
	return options
}

func normalizeOptions(options *Options) (Options, error) {
	normalized := Options{}
	if options != nil {
		normalized = *options
	}
	if normalized.MaxIterations < 0 {
		return Options{}, fmt.Errorf("nonlinear: %w: max iterations must be non-negative", adjust.ErrInvalidProblem)
	}
	if normalized.StepTolerance < 0 || math.IsNaN(normalized.StepTolerance) || math.IsInf(normalized.StepTolerance, 0) {
		return Options{}, fmt.Errorf("nonlinear: %w: step tolerance must be finite and non-negative", adjust.ErrInvalidProblem)
	}
	if normalized.ObjectiveTolerance < 0 || math.IsNaN(normalized.ObjectiveTolerance) || math.IsInf(normalized.ObjectiveTolerance, 0) {
		return Options{}, fmt.Errorf("nonlinear: %w: objective tolerance must be finite and non-negative", adjust.ErrInvalidProblem)
	}
	return normalized.withDefaults(), nil
}

// Result contains the converged state and the last linear adjustment.
type Result struct {
	State      []float64
	Adjustment *adjust.Result
	Iterations int
	Converged  bool
}

// Solve performs Gauss-Newton iteration from initial.
func Solve(initial []float64, model Model, options *Options) (*Result, error) {
	return SolveContext(context.Background(), initial, model, options)
}

// SolveContext performs Gauss-Newton iteration with cooperative cancellation.
// A running Model.Linearize call cannot be interrupted by this package, but
// cancellation is observed before and after every linearization and throughout
// the underlying batch solve.
func SolveContext(ctx context.Context, initial []float64, model Model, options *Options) (*Result, error) {
	if ctx == nil {
		return nil, fmt.Errorf("nonlinear: nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	if len(initial) == 0 {
		return nil, fmt.Errorf("nonlinear: %w: initial state is empty", adjust.ErrInvalidProblem)
	}
	if nilInterface(model) {
		return nil, fmt.Errorf("nonlinear: %w: model is nil", adjust.ErrInvalidProblem)
	}
	for i, value := range initial {
		if math.IsNaN(value) || math.IsInf(value, 0) {
			return nil, fmt.Errorf("nonlinear: %w: initial state element %d is not finite", adjust.ErrInvalidProblem, i)
		}
	}
	opts, err := normalizeOptions(options)
	if err != nil {
		return nil, err
	}
	state := append([]float64(nil), initial...)
	previousObjective := math.Inf(1)
	result := &Result{State: state}

	for iteration := 1; iteration <= opts.MaxIterations; iteration++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		problem, err := model.Linearize(append([]float64(nil), state...))
		if err != nil {
			return nil, fmt.Errorf("nonlinear: linearize iteration %d: %w", iteration, err)
		}
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		if problem.ParameterCount != len(state) {
			return nil, fmt.Errorf("nonlinear: %w: model returned %d parameters for a %d-element state", adjust.ErrInvalidProblem, problem.ParameterCount, len(state))
		}
		adjustment, err := batch.SolveContext(ctx, problem, opts.Batch)
		if err != nil {
			return nil, fmt.Errorf("nonlinear: iteration %d: %w", iteration, err)
		}
		maxStep := 0.0
		for i, correction := range adjustment.Delta {
			state[i] += correction
			maxStep = math.Max(maxStep, math.Abs(correction))
		}
		objectiveChange := math.Abs(previousObjective-adjustment.Objective) / math.Max(1, math.Abs(previousObjective))
		result.State = state
		result.Adjustment = adjustment
		result.Iterations = iteration
		objectiveStalled := iteration > 1 && objectiveChange <= opts.ObjectiveTolerance
		if maxStep <= opts.StepTolerance || (objectiveStalled && maxStep <= 100*opts.StepTolerance) {
			result.Converged = true
			return result, nil
		}
		previousObjective = adjustment.Objective
	}
	return result, nil
}

func nilInterface(value any) bool {
	if value == nil {
		return true
	}
	reflected := reflect.ValueOf(value)
	switch reflected.Kind() {
	case reflect.Chan, reflect.Func, reflect.Interface, reflect.Map, reflect.Pointer, reflect.Slice:
		return reflected.IsNil()
	default:
		return false
	}
}
