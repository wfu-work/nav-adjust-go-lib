package network

import (
	"errors"
	"fmt"
	"math"
	"sort"

	"github.com/wfu-work/nav-adjust-go-lib/batch"
	"gonum.org/v1/gonum/mat"
)

// ValidateENUNetwork validates identifiers, covariance matrices and network
// connectivity without running an adjustment.
func ValidateENUNetwork(problem ENUNetworkProblem, options *ENUNetworkOptions) error {
	opts, err := normalizeNetworkOptions(options)
	if err != nil {
		return err
	}
	_, err = compileENUNetwork(problem, opts)
	return err
}

// Validate checks identifiers, covariance matrices and network connectivity
// without running an adjustment.
func Validate(problem Problem, options *Options) error {
	return ValidateENUNetwork(problem, options)
}

func validateMatrix3(covariance Matrix3, options batch.Options, field, id string) error {
	symmetryTolerance := options.SymmetryTolerance
	if symmetryTolerance <= 0 {
		symmetryTolerance = 1e-12
	}
	minimumVariance := options.MinVariance
	if minimumVariance <= 0 {
		minimumVariance = 1e-20
	}
	matrix := mat.NewSymDense(3, nil)
	for row := range 3 {
		if !finite(covariance.At(row, row)) || covariance.At(row, row) < minimumVariance {
			return validationError(ErrInvalidCovariance, field, id, fmt.Sprintf("variance at component %d must be at least %g", row, minimumVariance))
		}
		for column := range 3 {
			value := covariance.At(row, column)
			if !finite(value) {
				return validationError(ErrInvalidCovariance, field, id, "contains a non-finite value")
			}
			if math.Abs(value-covariance.At(column, row)) > symmetryTolerance {
				return validationError(ErrInvalidCovariance, field, id, "must be symmetric")
			}
			if column <= row {
				matrix.SetSym(row, column, value)
			}
		}
	}
	var chol mat.Cholesky
	if !chol.Factorize(matrix) {
		return validationError(errors.Join(ErrInvalidCovariance, ErrNotPositiveDefinite), field, id, "must be positive definite")
	}
	condition := chol.Cond()
	if math.IsInf(condition, 1) || condition > 1e16 {
		return validationError(
			errors.Join(ErrInvalidCovariance, ErrNotPositiveDefinite), field, id,
			fmt.Sprintf("must be numerically positive definite; condition number %.6g is too large", condition),
		)
	}
	return nil
}

func datumComponents(
	order []Station,
	stations map[string]Station,
	adjacency map[string][]string,
	priorStations map[string]bool,
	allowFree bool,
) ([][]string, error) {
	visited := make(map[string]bool, len(stations))
	freeComponents := make([][]string, 0)
	for _, first := range order {
		if visited[first.ID] {
			continue
		}
		queue := []string{first.ID}
		visited[first.ID] = true
		anchored := false
		component := make([]string, 0)
		for len(queue) > 0 {
			id := queue[0]
			queue = queue[1:]
			component = append(component, id)
			anchored = anchored || stations[id].Fixed || priorStations[id]
			for _, neighbor := range adjacency[id] {
				if !visited[neighbor] {
					visited[neighbor] = true
					queue = append(queue, neighbor)
				}
			}
		}
		if !anchored {
			sort.Strings(component)
			if !allowFree {
				return nil, validationError(ErrDisconnectedNetwork, "stations", component[0], "component has no fixed station or position prior: "+fmt.Sprint(component))
			}
			if len(component) == 1 && len(adjacency[component[0]]) == 0 {
				return nil, validationError(ErrDisconnectedNetwork, "stations", component[0], "isolated station has no baseline or position prior")
			}
			freeComponents = append(freeComponents, component)
		}
	}
	return freeComponents, nil
}

func finiteENU(value ENU) bool {
	return finite(value.East) && finite(value.North) && finite(value.Up)
}

func finite(value float64) bool {
	return !math.IsNaN(value) && !math.IsInf(value, 0)
}
