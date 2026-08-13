// Package model defines the input and output contracts for ENU network
// adjustment without depending on a numerical solver.
package model

// ENU is a three-dimensional vector or position in metres.
type ENU struct {
	East  float64 `json:"east"`
	North float64 `json:"north"`
	Up    float64 `json:"up"`
}

// Matrix is a row-major dense matrix used by the public API. It deliberately
// avoids exposing the numerical backend used by the library.
type Matrix struct {
	Rows int       `json:"rows"`
	Cols int       `json:"cols"`
	Data []float64 `json:"data"`
}

// At returns the matrix value at zero-based row and column.
func (matrix Matrix) At(row, column int) float64 {
	return matrix.Data[row*matrix.Cols+column]
}

// Matrix3 is a row-major 3-by-3 matrix. For observation input it represents
// the covariance of [East, North, Up].
type Matrix3 struct {
	Data [9]float64 `json:"data"`
}

// DiagonalMatrix3 creates a diagonal 3-by-3 matrix from variances.
func DiagonalMatrix3(eastVariance, northVariance, upVariance float64) Matrix3 {
	return Matrix3{Data: [9]float64{
		eastVariance, 0, 0,
		0, northVariance, 0,
		0, 0, upVariance,
	}}
}

// Matrix3FromStdDev creates a diagonal covariance matrix from standard
// deviations.
func Matrix3FromStdDev(east, north, up float64) Matrix3 {
	return DiagonalMatrix3(east*east, north*north, up*up)
}

// At returns the matrix value at zero-based row and column.
func (matrix Matrix3) At(row, column int) float64 {
	return matrix.Data[row*3+column]
}
