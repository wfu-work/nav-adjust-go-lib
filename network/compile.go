package network

import (
	"fmt"
	"math"

	"github.com/wfu-work/nav-adjust-go-lib/core"
)

type compiledNetwork struct {
	problem             core.Problem
	stations            map[string]Station
	stationIndex        map[string][3]int
	parameterKeys       []string
	publicOptions       networkOptions
	dummyIndex          int
	freeDatumComponents [][]string
}

func compileENUNetwork(input ENUNetworkProblem, options networkOptions) (compiledNetwork, error) {
	if len(input.Stations) == 0 {
		return compiledNetwork{}, validationError(nil, "stations", "", "at least one station is required")
	}
	if len(input.Baselines) == 0 && len(input.Priors) == 0 {
		return compiledNetwork{}, validationError(nil, "baselines", "", "at least one baseline or position prior is required")
	}
	if options.variance != nil && len(input.Baselines) == 0 {
		return compiledNetwork{}, validationError(nil, "options.variance_components", "", "requires at least one baseline")
	}
	compiled := compiledNetwork{
		stations:      make(map[string]Station, len(input.Stations)),
		stationIndex:  make(map[string][3]int, len(input.Stations)),
		publicOptions: options,
	}
	freeStations := make([]Station, 0, len(input.Stations))
	for _, station := range input.Stations {
		if station.ID == "" {
			return compiledNetwork{}, validationError(nil, "stations.id", "", "must not be empty")
		}
		if _, exists := compiled.stations[station.ID]; exists {
			return compiledNetwork{}, validationError(ErrDuplicateStation, "stations.id", station.ID, "must be unique")
		}
		if station.Fixed {
			if station.KnownENU == nil {
				return compiledNetwork{}, validationError(nil, "stations.known_enu", station.ID, "is required for a fixed station")
			}
			if !finiteENU(*station.KnownENU) {
				return compiledNetwork{}, validationError(nil, "stations.known_enu", station.ID, "contains a non-finite coordinate")
			}
		} else {
			if station.KnownENU != nil {
				return compiledNetwork{}, validationError(nil, "stations.known_enu", station.ID, "is only valid for fixed stations")
			}
			freeStations = append(freeStations, station)
		}
		compiled.stations[station.ID] = station
	}

	for _, station := range freeStations {
		start := len(compiled.parameterKeys)
		compiled.stationIndex[station.ID] = [3]int{start, start + 1, start + 2}
		compiled.parameterKeys = append(compiled.parameterKeys,
			"station:"+station.ID+":east",
			"station:"+station.ID+":north",
			"station:"+station.ID+":up",
		)
	}
	solverParameterKeys := append([]string(nil), compiled.parameterKeys...)
	needsDummy := len(freeStations) == 0
	if !needsDummy {
		for _, baseline := range input.Baselines {
			from, fromExists := compiled.stations[baseline.From]
			to, toExists := compiled.stations[baseline.To]
			if fromExists && toExists && from.Fixed && to.Fixed {
				needsDummy = true
				break
			}
		}
	}
	if needsDummy {
		// The numerical core requires a non-zero design row. A constrained dummy
		// lets fixed-to-fixed baselines participate in residual diagnostics.
		compiled.dummyIndex = len(solverParameterKeys)
		solverParameterKeys = append(solverParameterKeys, "internal:fixed-network")
	}
	compiled.problem = core.NewProblem(solverParameterKeys...)
	if needsDummy {
		compiled.problem.AddConstraint(core.FixCorrection(compiled.dummyIndex, 0))
	}

	baselineIDs := make(map[string]struct{}, len(input.Baselines))
	adjacency := make(map[string][]string, len(input.Stations))
	for id := range compiled.stations {
		adjacency[id] = nil
	}
	for baselineIndex, baseline := range input.Baselines {
		if baseline.ID == "" {
			return compiledNetwork{}, validationError(nil, "baselines.id", fmt.Sprintf("index:%d", baselineIndex), "must not be empty")
		}
		if _, exists := baselineIDs[baseline.ID]; exists {
			return compiledNetwork{}, validationError(ErrDuplicateBaseline, "baselines.id", baseline.ID, "must be unique")
		}
		baselineIDs[baseline.ID] = struct{}{}
		from, fromExists := compiled.stations[baseline.From]
		to, toExists := compiled.stations[baseline.To]
		if !fromExists {
			return compiledNetwork{}, validationError(ErrUnknownStation, "baselines.from", baseline.ID, fmt.Sprintf("station %q does not exist", baseline.From))
		}
		if !toExists {
			return compiledNetwork{}, validationError(ErrUnknownStation, "baselines.to", baseline.ID, fmt.Sprintf("station %q does not exist", baseline.To))
		}
		if baseline.From == baseline.To {
			return compiledNetwork{}, validationError(nil, "baselines.to", baseline.ID, "must differ from from")
		}
		if !finiteENU(baseline.Vector) {
			return compiledNetwork{}, validationError(nil, "baselines.vector", baseline.ID, "contains a non-finite value")
		}
		if err := validateMatrix3(baseline.Covariance, options.batch, "baselines.covariance", baseline.ID); err != nil {
			return compiledNetwork{}, err
		}
		adjacency[baseline.From] = append(adjacency[baseline.From], baseline.To)
		adjacency[baseline.To] = append(adjacency[baseline.To], baseline.From)

		fromPosition := stationPosition(from)
		toPosition := stationPosition(to)
		observed := [3]float64{baseline.Vector.East, baseline.Vector.North, baseline.Vector.Up}
		fromKnown := [3]float64{fromPosition.East, fromPosition.North, fromPosition.Up}
		toKnown := [3]float64{toPosition.East, toPosition.North, toPosition.Up}
		rowIndexes := [3]int{}
		components := [3]string{"east", "north", "up"}
		for component := range 3 {
			terms := make([]core.Term, 0, 2)
			if from.Fixed && to.Fixed {
				terms = append(terms, core.T(compiled.dummyIndex, 1))
			}
			if !from.Fixed {
				terms = append(terms, core.T(compiled.stationIndex[from.ID][component], -1))
			}
			if !to.Fixed {
				terms = append(terms, core.T(compiled.stationIndex[to.ID][component], 1))
			}
			rowIndexes[component] = compiled.problem.AddEquation(core.Equation{
				ID: baseline.ID + ":" + components[component], Group: baseline.Group,
				Terms: terms, Misclosure: observed[component] + fromKnown[component] - toKnown[component],
			})
		}
		compiled.problem.CovarianceBlocks = append(compiled.problem.CovarianceBlocks, core.CovarianceBlock{
			ID: baseline.ID, RowIndexes: rowIndexes[:],
			Covariance: append([]float64(nil), baseline.Covariance.Data[:]...),
		})
	}

	priorStations, err := compiled.addPositionPriors(input.Priors, options)
	if err != nil {
		return compiledNetwork{}, err
	}
	freeComponents, err := datumComponents(
		input.Stations, compiled.stations, adjacency, priorStations, options.datum == DatumFreeCentroid,
	)
	if err != nil {
		return compiledNetwork{}, err
	}
	compiled.freeDatumComponents = freeComponents
	compiled.addFreeDatumConstraints()
	return compiled, nil
}

func (compiled *compiledNetwork) addFreeDatumConstraints() {
	components := [3]string{"east", "north", "up"}
	for componentIndex, stationIDs := range compiled.freeDatumComponents {
		coefficient := 1 / math.Sqrt(float64(len(stationIDs)))
		for axis := range 3 {
			terms := make([]core.Term, 0, len(stationIDs))
			for _, stationID := range stationIDs {
				terms = append(terms, core.T(compiled.stationIndex[stationID][axis], coefficient))
			}
			compiled.problem.AddConstraint(core.ExactConstraint(
				fmt.Sprintf("datum:free-centroid:%d:%s", componentIndex, components[axis]), 0, terms...,
			))
		}
	}
}

func stationPosition(station Station) ENU {
	if station.KnownENU != nil {
		return *station.KnownENU
	}
	return ENU{}
}
