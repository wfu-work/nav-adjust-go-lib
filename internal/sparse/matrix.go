// Package sparse provides the small sparse linear-algebra surface needed by
// the adjustment solvers. It is internal so the storage format can evolve
// without becoming part of the public API.
package sparse

import (
	"fmt"
	"math"
	"sort"
)

// Builder incrementally assembles a symmetric matrix.
type Builder struct {
	size int
	rows []map[int]float64
}

// NewBuilder creates a symmetric matrix builder of the requested size.
func NewBuilder(size int) *Builder {
	rows := make([]map[int]float64, size)
	for i := range rows {
		rows[i] = make(map[int]float64)
	}
	return &Builder{size: size, rows: rows}
}

// AddSym adds value to both symmetric positions of the matrix.
func (builder *Builder) AddSym(row, column int, value float64) {
	if value == 0 {
		return
	}
	builder.rows[row][column] += value
	if row != column {
		builder.rows[column][row] += value
	}
}

// Matrix is an immutable compressed sparse row matrix.
type Matrix struct {
	size     int
	rowStart []int
	columns  []int
	values   []float64
	diagonal []float64
}

// Build finalizes the matrix and releases no references to the builder maps.
func (builder *Builder) Build() (*Matrix, error) {
	if builder == nil || builder.size <= 0 || len(builder.rows) != builder.size {
		return nil, fmt.Errorf("sparse: invalid builder")
	}
	result := &Matrix{
		size:     builder.size,
		rowStart: make([]int, builder.size+1),
		diagonal: make([]float64, builder.size),
	}
	for row, entries := range builder.rows {
		columns := make([]int, 0, len(entries))
		for column, value := range entries {
			if value != 0 {
				columns = append(columns, column)
			}
		}
		sort.Ints(columns)
		for _, column := range columns {
			value := entries[column]
			if math.IsNaN(value) || math.IsInf(value, 0) {
				return nil, fmt.Errorf("sparse: matrix contains a non-finite value")
			}
			result.columns = append(result.columns, column)
			result.values = append(result.values, value)
			if row == column {
				result.diagonal[row] = value
			}
		}
		result.rowStart[row+1] = len(result.columns)
	}
	return result, nil
}

// Size returns the matrix dimension.
func (matrix *Matrix) Size() int { return matrix.size }

// NNZ returns the number of explicitly stored entries. Both triangles are
// counted because multiplication uses a conventional CSR representation.
func (matrix *Matrix) NNZ() int { return len(matrix.values) }

// Diagonal returns a copy of the matrix diagonal.
func (matrix *Matrix) Diagonal() []float64 {
	return append([]float64(nil), matrix.diagonal...)
}

// MulVec stores matrix*x in destination. Destination and x must not alias.
func (matrix *Matrix) MulVec(destination, x []float64) {
	for row := 0; row < matrix.size; row++ {
		value := 0.0
		for offset := matrix.rowStart[row]; offset < matrix.rowStart[row+1]; offset++ {
			value += matrix.values[offset] * x[matrix.columns[offset]]
		}
		destination[row] = value
	}
}

// At returns one matrix entry, including an implicit zero.
func (matrix *Matrix) At(row, column int) float64 {
	start, end := matrix.rowStart[row], matrix.rowStart[row+1]
	offset := sort.Search(end-start, func(i int) bool {
		return matrix.columns[start+i] >= column
	})
	if start+offset < end && matrix.columns[start+offset] == column {
		return matrix.values[start+offset]
	}
	return 0
}
