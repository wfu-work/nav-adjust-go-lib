// Command spp demonstrates an iterative ECEF single-point-position adjustment
// with a receiver clock bias expressed in metres.
package main

import (
	"fmt"
	"math"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
	"github.com/wfu-work/nav-adjust-go-lib/nonlinear"
)

type satellite struct {
	position    [3]float64
	pseudorange float64
	variance    float64
}

func main() {
	trueReceiver := [4]float64{-2267800, 5009300, 3221000, 75000}
	positions := [][3]float64{
		{15600e3, 7540e3, 20140e3},
		{18760e3, 2750e3, 18610e3},
		{17610e3, 14630e3, 13480e3},
		{19170e3, 610e3, 18390e3},
		{17800e3, -15200e3, 15100e3},
		{-13400e3, 21300e3, 9000e3},
	}
	noise := []float64{1.2, -0.8, 0.4, -1.0, 0.7, -0.2}
	satellites := make([]satellite, len(positions))
	for i, position := range positions {
		satellites[i] = satellite{
			position:    position,
			pseudorange: distance(trueReceiver[:3], position[:]) + trueReceiver[3] + noise[i],
			variance:    4,
		}
	}

	model := nonlinear.ModelFunc(func(state []float64) (adjust.Problem, error) {
		problem := adjust.NewProblem("ecef-x", "ecef-y", "ecef-z", "clock-metres")
		for i, sat := range satellites {
			dx := state[0] - sat.position[0]
			dy := state[1] - sat.position[1]
			dz := state[2] - sat.position[2]
			rho := math.Sqrt(dx*dx + dy*dy + dz*dz)
			problem.AddEquation(adjust.Equation{
				ID: fmt.Sprintf("G%02d-C1C", i+1), Group: "GPS-C1C",
				Terms: []adjust.Term{
					adjust.T(0, dx/rho), adjust.T(1, dy/rho),
					adjust.T(2, dz/rho), adjust.T(3, 1),
				},
				Misclosure: sat.pseudorange - (rho + state[3]),
				Variance:   sat.variance,
			})
		}
		return problem, nil
	})

	initial := []float64{-2266800, 5008300, 3221500, 0}
	result, err := nonlinear.Solve(initial, model, nil)
	if err != nil {
		panic(err)
	}
	fmt.Printf("converged=%v iterations=%d\n", result.Converged, result.Iterations)
	fmt.Printf("ECEF = %.3f %.3f %.3f m\n", result.State[0], result.State[1], result.State[2])
	fmt.Printf("clock = %.3f m, sigma0 = %.3f\n", result.State[3], result.Adjustment.Sigma0)
}

func distance(a, b []float64) float64 {
	dx := a[0] - b[0]
	dy := a[1] - b[1]
	dz := a[2] - b[2]
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}
