package network

import (
	"context"
	"fmt"
	"math"
	"sort"

	"github.com/wfu-work/nav-adjust-go-lib/batch"
	"github.com/wfu-work/nav-adjust-go-lib/core"
	"gonum.org/v1/gonum/mat"
)

type covarianceReader interface {
	At(row, column int) float64
}

type selectedCovariance map[[2]int]float64

func (covariance selectedCovariance) At(row, column int) float64 {
	if row > column {
		row, column = column, row
	}
	return covariance[[2]int{row, column}]
}

func queryNetworkCovariance(ctx context.Context, solution *batch.Solution, problem core.Problem, stationIndexes map[string][3]int) (selectedCovariance, error) {
	pairSet := make(map[[2]int]struct{})
	for _, indexes := range stationIndexes {
		addIndexPairs(pairSet, indexes[:])
	}
	for _, block := range problem.CovarianceBlocks {
		parameters := make(map[int]struct{})
		for _, rowIndex := range block.RowIndexes {
			for _, term := range problem.Equations[rowIndex].Terms {
				parameters[term.Parameter] = struct{}{}
			}
		}
		indexes := make([]int, 0, len(parameters))
		for parameter := range parameters {
			indexes = append(indexes, parameter)
		}
		sort.Ints(indexes)
		addIndexPairs(pairSet, indexes)
	}
	pairs := make([][2]int, 0, len(pairSet))
	for pair := range pairSet {
		pairs = append(pairs, pair)
	}
	sort.Slice(pairs, func(i, j int) bool {
		if pairs[i][1] != pairs[j][1] {
			return pairs[i][1] < pairs[j][1]
		}
		return pairs[i][0] < pairs[j][0]
	})
	values, err := solution.FormalCovarianceValuesContext(ctx, pairs)
	if err != nil {
		return nil, fmt.Errorf("adjust: query local covariance: %w", err)
	}
	result := make(selectedCovariance, len(pairs))
	for i, pair := range pairs {
		result[pair] = values[i]
	}
	return result, nil
}

func addIndexPairs(destination map[[2]int]struct{}, indexes []int) {
	for r, row := range indexes {
		for c := 0; c <= r; c++ {
			column := indexes[c]
			if row > column {
				destination[[2]int{column, row}] = struct{}{}
			} else {
				destination[[2]int{row, column}] = struct{}{}
			}
		}
	}
}

func baselineResidualCovariance(problem core.Problem, covariance covarianceReader, baselineIndex int) Matrix3 {
	return observationResidualCovariance(problem, covariance, baselineIndex)
}

func observationResidualCovariance(problem core.Problem, covariance covarianceReader, blockIndex int) Matrix3 {
	block := problem.CovarianceBlocks[blockIndex]
	result := Matrix3{}
	for row := range 3 {
		left := problem.Equations[block.RowIndexes[row]].Terms
		for column := 0; column <= row; column++ {
			right := problem.Equations[block.RowIndexes[column]].Terms
			value := block.Covariance[row*3+column] - sparseCrossCovariance(left, covariance, right)
			result.Data[row*3+column] = value
			result.Data[column*3+row] = value
		}
	}
	return result
}

func sparseCrossCovariance(left []core.Term, covariance covarianceReader, right []core.Term) float64 {
	value := 0.0
	for _, leftTerm := range left {
		for _, rightTerm := range right {
			value += leftTerm.Coefficient * covariance.At(leftTerm.Parameter, rightTerm.Parameter) * rightTerm.Coefficient
		}
	}
	return value
}

func publicMatrix(matrix mat.Matrix, size int) Matrix {
	result := Matrix{Rows: size, Cols: size, Data: make([]float64, size*size)}
	for row := range size {
		for column := range size {
			result.Data[row*size+column] = matrix.At(row, column)
		}
	}
	return result
}

func covarianceBlock3(covariance covarianceReader, indexes [3]int) Matrix3 {
	result := Matrix3{}
	for row := range 3 {
		for column := range 3 {
			result.Data[row*3+column] = covariance.At(indexes[row], indexes[column])
		}
	}
	return result
}

func observationDiagnostics(formalResidual Matrix3, observationCovariance []float64, residual ENU, sigma0 float64) (ENU, ENU, error) {
	formal := mat.NewSymDense(3, nil)
	observation := mat.NewSymDense(3, nil)
	for row := range 3 {
		for column := 0; column <= row; column++ {
			formal.SetSym(row, column, formalResidual.At(row, column))
			observation.SetSym(row, column, observationCovariance[row*3+column])
		}
	}
	var factor mat.Cholesky
	if !factor.Factorize(observation) {
		return ENU{}, ENU{}, fmt.Errorf("%w: observation covariance", ErrNotPositiveDefinite)
	}
	var weighted mat.Dense
	if err := factor.SolveTo(&weighted, formal); err != nil {
		return ENU{}, ENU{}, fmt.Errorf("adjust: residual redundancy: %w", err)
	}
	components := [3]float64{residual.East, residual.North, residual.Up}
	standardized := [3]float64{}
	redundancy := [3]float64{}
	for component := range 3 {
		variance := nonnegativeCovariance(formal.At(component, component), observation.At(component, component))
		if variance > 0 && sigma0 > 0 {
			standardized[component] = components[component] / (sigma0 * math.Sqrt(variance))
		}
		redundancy[component] = math.Max(0, math.Min(1, weighted.At(component, component)))
	}
	return ENU{East: standardized[0], North: standardized[1], Up: standardized[2]},
		ENU{East: redundancy[0], North: redundancy[1], Up: redundancy[2]}, nil
}

func nonnegativeCovariance(value, scale float64) float64 {
	if value < 0 && value > -1e-10*math.Max(1, scale) {
		return 0
	}
	return value
}

func scaleMatrix3(matrix Matrix3, scale float64) Matrix3 {
	for i := range matrix.Data {
		matrix.Data[i] *= scale
	}
	return matrix
}
