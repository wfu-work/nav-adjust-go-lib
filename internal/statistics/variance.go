package statistics

import (
	"errors"
	"fmt"
	"math"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
)

var (
	// ErrInvalidVarianceComponent indicates non-finite or inconsistent
	// component statistics.
	ErrInvalidVarianceComponent = errors.New("adjust: invalid variance component")
)

// VarianceComponent contains the statistics required to update one covariance
// scale.
type VarianceComponent struct {
	ID               string
	Scale            float64
	Objective        float64
	Redundancy       float64
	BaselineCount    int
	ObservationCount int
}

// VarianceOptions controls one group-scale update.
type VarianceOptions struct {
	Tolerance         float64
	MinScale          float64
	MaxScale          float64
	MinimumRedundancy float64
}

// VarianceUpdate is the result of one variance-component update.
type VarianceUpdate struct {
	Components        []VarianceComponent
	MaxRelativeChange float64
	Converged         bool
}

// UpdateGroupScales applies one group-scale update without mutating components.
func UpdateGroupScales(components []VarianceComponent, options VarianceOptions) (VarianceUpdate, error) {
	if err := validateVarianceOptions(options); err != nil {
		return VarianceUpdate{}, err
	}
	result := VarianceUpdate{Components: make([]VarianceComponent, len(components)), Converged: true}
	for i, component := range components {
		if !positiveFinite(component.Scale) {
			return VarianceUpdate{}, fmt.Errorf("%w %q: scale must be positive and finite", ErrInvalidVarianceComponent, component.ID)
		}
		if component.Objective < 0 || !finite(component.Objective) {
			return VarianceUpdate{}, fmt.Errorf("%w %q: objective must be non-negative and finite", ErrInvalidVarianceComponent, component.ID)
		}
		if !finite(component.Redundancy) || component.Redundancy < options.MinimumRedundancy {
			return VarianceUpdate{}, fmt.Errorf("%w for group %q: got %.6g, need at least %.6g",
				adjust.ErrInsufficientRedundancy, component.ID, component.Redundancy, options.MinimumRedundancy)
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
			return VarianceUpdate{}, fmt.Errorf("%w %q: updated scale is not positive and finite", ErrInvalidVarianceComponent, component.ID)
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

func validateVarianceOptions(options VarianceOptions) error {
	if options.Tolerance <= 0 || !finite(options.Tolerance) {
		return fmt.Errorf("%w: tolerance must be positive and finite", ErrInvalidVarianceComponent)
	}
	if !positiveFinite(options.MinScale) {
		return fmt.Errorf("%w: minimum scale must be positive and finite", ErrInvalidVarianceComponent)
	}
	if !positiveFinite(options.MaxScale) || options.MaxScale < options.MinScale {
		return fmt.Errorf("%w: maximum scale must be finite and no smaller than minimum scale", ErrInvalidVarianceComponent)
	}
	if options.MinimumRedundancy <= 0 || !finite(options.MinimumRedundancy) {
		return fmt.Errorf("%w: minimum redundancy must be positive and finite", ErrInvalidVarianceComponent)
	}
	return nil
}

func positiveFinite(value float64) bool {
	return value > 0 && finite(value)
}

func finite(value float64) bool {
	return !math.IsNaN(value) && !math.IsInf(value, 0)
}
