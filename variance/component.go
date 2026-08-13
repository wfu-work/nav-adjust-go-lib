// Package variance provides numerical helpers for variance-component
// estimation. It is independent of the ENU network data model.
package variance

import (
	"errors"
	"fmt"
	"math"
)

var (
	// ErrInvalidComponent indicates non-finite or otherwise inconsistent
	// component statistics.
	ErrInvalidComponent = errors.New("adjust: invalid variance component")
	// ErrInsufficientRedundancy indicates that a component cannot be estimated
	// independently from the supplied observations.
	ErrInsufficientRedundancy = errors.New("adjust: insufficient variance-component redundancy")
)

// Component contains the statistics required to update one covariance scale.
// Scale multiplies the component's input covariance matrix. Objective is its
// weighted residual sum of squares and Redundancy is its effective degrees of
// freedom under the current stochastic model.
type Component struct {
	ID               string
	Scale            float64
	Objective        float64
	Redundancy       float64
	BaselineCount    int
	ObservationCount int
}

// Options controls one group-scale update.
type Options struct {
	Tolerance         float64
	MinScale          float64
	MaxScale          float64
	MinimumRedundancy float64
}

// Update is the result of one variance-component update.
type Update struct {
	Components        []Component
	MaxRelativeChange float64
	Converged         bool
}

// UpdateGroupScales applies an iterative group-scale update:
//
//	new scale = old scale * objective / redundancy
//
// The result is clamped to the configured scale interval. The function does
// not mutate components.
func UpdateGroupScales(components []Component, options Options) (Update, error) {
	if err := validateOptions(options); err != nil {
		return Update{}, err
	}
	result := Update{Components: make([]Component, len(components)), Converged: true}
	for i, component := range components {
		if !positiveFinite(component.Scale) {
			return Update{}, fmt.Errorf("%w %q: scale must be positive and finite", ErrInvalidComponent, component.ID)
		}
		if component.Objective < 0 || !finite(component.Objective) {
			return Update{}, fmt.Errorf("%w %q: objective must be non-negative and finite", ErrInvalidComponent, component.ID)
		}
		if !finite(component.Redundancy) || component.Redundancy < options.MinimumRedundancy {
			return Update{}, fmt.Errorf("%w for group %q: got %.6g, need at least %.6g",
				ErrInsufficientRedundancy, component.ID, component.Redundancy, options.MinimumRedundancy)
		}
		factor := component.Objective / component.Redundancy
		nextScale := component.Scale * factor
		if nextScale < options.MinScale {
			nextScale = options.MinScale
		}
		if nextScale > options.MaxScale {
			nextScale = options.MaxScale
		}
		if !positiveFinite(nextScale) {
			return Update{}, fmt.Errorf("%w %q: updated scale is not positive and finite", ErrInvalidComponent, component.ID)
		}
		change := math.Abs(nextScale/component.Scale - 1)
		result.MaxRelativeChange = math.Max(result.MaxRelativeChange, change)
		if change > options.Tolerance {
			result.Converged = false
		}
		component.Scale = nextScale
		result.Components[i] = component
	}
	return result, nil
}

func validateOptions(options Options) error {
	if options.Tolerance <= 0 || !finite(options.Tolerance) {
		return fmt.Errorf("%w: tolerance must be positive and finite", ErrInvalidComponent)
	}
	if !positiveFinite(options.MinScale) {
		return fmt.Errorf("%w: minimum scale must be positive and finite", ErrInvalidComponent)
	}
	if !positiveFinite(options.MaxScale) || options.MaxScale < options.MinScale {
		return fmt.Errorf("%w: maximum scale must be finite and no smaller than minimum scale", ErrInvalidComponent)
	}
	if options.MinimumRedundancy <= 0 || !finite(options.MinimumRedundancy) {
		return fmt.Errorf("%w: minimum redundancy must be positive and finite", ErrInvalidComponent)
	}
	return nil
}

func positiveFinite(value float64) bool {
	return value > 0 && finite(value)
}

func finite(value float64) bool {
	return !math.IsNaN(value) && !math.IsInf(value, 0)
}
