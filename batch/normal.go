package batch

import (
	"context"
	"fmt"

	"github.com/wfu-work/nav-adjust-go-lib/internal/sparse"
	"gonum.org/v1/gonum/mat"
)

type normalEquation struct {
	n *mat.SymDense
	u *mat.VecDense
}

type sparseNormalEquation struct {
	n *sparse.Matrix
	u []float64
}

type normalAccumulator interface {
	addSym(row, column int, value float64)
	addRHS(index int, value float64)
}

type denseAccumulator struct{ normal normalEquation }

func (accumulator denseAccumulator) addSym(row, column int, value float64) {
	accumulator.normal.n.SetSym(row, column, accumulator.normal.n.At(row, column)+value)
}

func (accumulator denseAccumulator) addRHS(index int, value float64) {
	accumulator.normal.u.SetVec(index, accumulator.normal.u.AtVec(index)+value)
}

type sparseAccumulator struct {
	builder *sparse.Builder
	u       []float64
}

func (accumulator sparseAccumulator) addSym(row, column int, value float64) {
	accumulator.builder.AddSym(row, column, value)
}

func (accumulator sparseAccumulator) addRHS(index int, value float64) {
	accumulator.u[index] += value
}

func buildNormal(ctx context.Context, compiled compiledProblem) (normalEquation, error) {
	nx := compiled.parameters
	normal := normalEquation{
		n: mat.NewSymDense(nx, nil),
		u: mat.NewVecDense(nx, nil),
	}
	accumulator := denseAccumulator{normal: normal}
	if err := accumulateNormal(ctx, accumulator, compiled); err != nil {
		return normalEquation{}, err
	}
	return normal, nil
}

func buildSparseNormal(ctx context.Context, compiled compiledProblem) (sparseNormalEquation, error) {
	accumulator := sparseAccumulator{
		builder: sparse.NewBuilder(compiled.parameters),
		u:       make([]float64, compiled.parameters),
	}
	if err := accumulateNormal(ctx, accumulator, compiled); err != nil {
		return sparseNormalEquation{}, err
	}
	matrix, err := accumulator.builder.Build()
	if err != nil {
		return sparseNormalEquation{}, fmt.Errorf("batch: build sparse normal equation: %w", err)
	}
	return sparseNormalEquation{n: matrix, u: accumulator.u}, nil
}

func accumulateNormal(ctx context.Context, accumulator normalAccumulator, compiled compiledProblem) error {
	for blockIndex := range compiled.blocks {
		if err := contextError(ctx); err != nil {
			return err
		}
		if err := accumulateBlock(accumulator, compiled.rows, &compiled.blocks[blockIndex]); err != nil {
			return err
		}
	}
	for i, row := range compiled.rows {
		if err := contextError(ctx); err != nil {
			return err
		}
		if compiled.blocked[i] {
			continue
		}
		accumulateIndependent(accumulator, row)
	}
	return nil
}

func accumulateIndependent(normal normalAccumulator, row stochasticRow) {
	weight := 1 / row.variance
	for column, columnTerm := range row.terms {
		c := columnTerm.Parameter
		hc := columnTerm.Coefficient
		normal.addRHS(c, hc*weight*row.w)
		for rowIndex := 0; rowIndex <= column; rowIndex++ {
			rowTerm := row.terms[rowIndex]
			r := rowTerm.Parameter
			normal.addSym(r, c, rowTerm.Coefficient*weight*hc)
		}
	}
}

func accumulateBlock(normal normalAccumulator, rows []stochasticRow, block *stochasticBlock) error {
	n := len(block.rowIndexes)
	parameterCount := len(block.parameters)
	parameterColumns := make(map[int]int, parameterCount)
	for column, parameter := range block.parameters {
		parameterColumns[parameter] = column
	}
	h := mat.NewDense(n, parameterCount, nil)
	w := mat.NewDense(n, 1, nil)
	for r, rowIndex := range block.rowIndexes {
		for _, term := range rows[rowIndex].terms {
			h.Set(r, parameterColumns[term.Parameter], term.Coefficient)
		}
		w.Set(r, 0, rows[rowIndex].w)
	}
	weightedH := mat.NewDense(n, parameterCount, nil)
	if err := block.cholesky.SolveTo(weightedH, h); err != nil {
		return fmt.Errorf("batch: solve covariance block %q: %w", block.id, err)
	}
	weightedW := mat.NewDense(n, 1, nil)
	if err := block.cholesky.SolveTo(weightedW, w); err != nil {
		return fmt.Errorf("batch: solve covariance block %q misclosure: %w", block.id, err)
	}
	for c, globalC := range block.parameters {
		for i := 0; i < n; i++ {
			normal.addRHS(globalC, h.At(i, c)*weightedW.At(i, 0))
		}
		for r := 0; r <= c; r++ {
			globalR := block.parameters[r]
			value := 0.0
			for i := 0; i < n; i++ {
				value += h.At(i, r) * weightedH.At(i, c)
			}
			normal.addSym(globalR, globalC, value)
		}
	}
	return nil
}
