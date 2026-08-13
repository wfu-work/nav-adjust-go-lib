package network

import (
	"fmt"

	"github.com/wfu-work/nav-adjust-go-lib/core"
)

// addPositionPriors appends stochastic control coordinates after all baseline
// observations. Keeping this ordering stable lets robust weighting operate on
// baseline blocks without changing the control-coordinate random model.
func (compiled *compiledNetwork) addPositionPriors(priors []PositionPrior, options networkOptions) (map[string]bool, error) {
	priorIDs := make(map[string]struct{}, len(priors))
	priorStations := make(map[string]bool, len(priors))
	for priorIndex, prior := range priors {
		if prior.ID == "" {
			return nil, validationError(nil, "priors.id", fmt.Sprintf("index:%d", priorIndex), "must not be empty")
		}
		if _, exists := priorIDs[prior.ID]; exists {
			return nil, validationError(ErrDuplicatePrior, "priors.id", prior.ID, "must be unique")
		}
		priorIDs[prior.ID] = struct{}{}

		station, exists := compiled.stations[prior.StationID]
		if !exists {
			return nil, validationError(ErrUnknownStation, "priors.station_id", prior.ID, fmt.Sprintf("station %q does not exist", prior.StationID))
		}
		if station.Fixed {
			return nil, validationError(nil, "priors.station_id", prior.ID, "must reference a free station; fixed stations are exact")
		}
		if !finiteENU(prior.Position) {
			return nil, validationError(nil, "priors.position", prior.ID, "contains a non-finite coordinate")
		}
		if err := validateMatrix3(prior.Covariance, options.batch, "priors.covariance", prior.ID); err != nil {
			return nil, err
		}

		indexes := compiled.stationIndex[prior.StationID]
		observed := [3]float64{prior.Position.East, prior.Position.North, prior.Position.Up}
		components := [3]string{"east", "north", "up"}
		rowIndexes := [3]int{}
		for component := range 3 {
			rowIndexes[component] = compiled.problem.AddEquation(core.Equation{
				ID:         "prior:" + prior.ID + ":" + components[component],
				Group:      "position-prior",
				Terms:      []core.Term{core.T(indexes[component], 1)},
				Misclosure: observed[component],
			})
		}
		compiled.problem.CovarianceBlocks = append(compiled.problem.CovarianceBlocks, core.CovarianceBlock{
			ID:         "prior:" + prior.ID,
			RowIndexes: rowIndexes[:],
			Covariance: append([]float64(nil), prior.Covariance.Data[:]...),
		})
		priorStations[prior.StationID] = true
	}
	return priorStations, nil
}
