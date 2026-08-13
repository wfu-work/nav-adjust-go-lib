package main

import (
	"fmt"

	"github.com/wfu-work/nav-adjust-go-lib/model"
	"github.com/wfu-work/nav-adjust-go-lib/network"
)

func main() {
	datum := model.ENU{}
	covariance := model.Matrix3FromStdDev(0.01, 0.01, 0.02)
	problem := model.Problem{
		Name: "demo-network",
		Stations: []model.Station{
			{ID: "A", Name: "datum", Fixed: true, KnownENU: &datum},
			{ID: "B"},
			{ID: "C"},
		},
		Baselines: []model.Baseline{
			{ID: "AB", From: "A", To: "B", Vector: model.ENU{East: 10.002, North: 0.001, Up: 0.003}, Covariance: covariance},
			{ID: "AC", From: "A", To: "C", Vector: model.ENU{East: 0.001, North: 9.998, Up: -0.002}, Covariance: covariance},
			{ID: "BC", From: "B", To: "C", Vector: model.ENU{East: -10.000, North: 10.001, Up: -0.004}, Covariance: covariance},
		},
	}

	if err := network.Validate(problem, nil); err != nil {
		panic(err)
	}
	result, err := network.Solve(problem, nil)
	if err != nil {
		panic(err)
	}
	for _, station := range result.Stations {
		fmt.Printf("%s: E=%.4f N=%.4f U=%.4f\n",
			station.ID, station.Position.East, station.Position.North, station.Position.Up)
	}
	fmt.Printf("sigma0=%.4f, dof=%d\n", result.Diagnostics.Sigma0, result.Diagnostics.DegreesOfFreedom)
}
