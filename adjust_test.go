package adjust_test

import (
	"testing"

	adjust "github.com/wfu-work/nav-adjust-go-lib"
)

func TestRootFacadeSolvesNetwork(t *testing.T) {
	origin := adjust.ENU{}
	problem := adjust.ENUNetworkProblem{
		Stations: []adjust.Station{
			{ID: "A", Fixed: true, KnownENU: &origin},
			{ID: "B"},
		},
		Baselines: []adjust.ENUBaseline{
			{
				ID: "AB", From: "A", To: "B",
				Vector:     adjust.ENU{East: 1, North: 2, Up: 3},
				Covariance: adjust.Matrix3FromStdDev(0.01, 0.01, 0.02),
			},
		},
	}
	if err := adjust.ValidateENUNetwork(problem, nil); err != nil {
		t.Fatal(err)
	}
	result, err := adjust.SolveENUNetwork(problem, nil)
	if err != nil {
		t.Fatal(err)
	}
	if got := result.Stations[1].Position; got != (adjust.ENU{East: 1, North: 2, Up: 3}) {
		t.Fatalf("position = %+v", got)
	}
}

func TestRootFacadeExposesPositionPrior(t *testing.T) {
	problem := adjust.ENUNetworkProblem{
		Stations: []adjust.Station{{ID: "A"}},
		Priors: []adjust.PositionPrior{{
			ID: "control-a", StationID: "A",
			Position:   adjust.ENU{East: 1, North: 2, Up: 3},
			Covariance: adjust.DiagonalMatrix3(1, 1, 1),
		}},
	}
	result, err := adjust.SolveENUNetwork(problem, nil)
	if err != nil {
		t.Fatal(err)
	}
	if len(result.Priors) != 1 || result.Stations[0].Position != problem.Priors[0].Position {
		t.Fatalf("unexpected prior result: %+v", result)
	}
}

func TestRootFacadeExposesFreeCentroidDatum(t *testing.T) {
	problem := adjust.ENUNetworkProblem{
		Stations: []adjust.Station{{ID: "A"}, {ID: "B"}},
		Baselines: []adjust.ENUBaseline{{
			ID: "AB", From: "A", To: "B", Vector: adjust.ENU{East: 2},
			Covariance: adjust.DiagonalMatrix3(1, 1, 1),
		}},
	}
	result, err := adjust.SolveENUNetwork(problem, &adjust.ENUNetworkOptions{Datum: adjust.DatumFreeCentroid})
	if err != nil {
		t.Fatal(err)
	}
	if result.Stations[0].Position.East != -1 || result.Stations[1].Position.East != 1 ||
		result.Diagnostics.DatumMode != adjust.DatumFreeCentroid {
		t.Fatalf("unexpected free-datum facade result: %+v", result)
	}
}
