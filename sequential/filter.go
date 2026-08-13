// Package sequential implements linear Kalman prediction and measurement
// update for epoch-by-epoch GNSS estimators.
package sequential

import (
	"fmt"
	"math"
	"reflect"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
	"gonum.org/v1/gonum/mat"
)

// State is a Kalman state vector and covariance.
type State struct {
	X *mat.VecDense
	P *mat.SymDense
}

// NewState copies x and covariance into an independent state.
func NewState(x []float64, covariance *mat.SymDense) (*State, error) {
	if len(x) == 0 || covariance == nil || covariance.SymmetricDim() != len(x) {
		return nil, fmt.Errorf("sequential: %w: state and covariance dimensions do not match", adjust.ErrInvalidProblem)
	}
	if err := validateFiniteVector(mat.NewVecDense(len(x), x), "state"); err != nil {
		return nil, err
	}
	if err := validateCovariance(covariance, len(x), "state covariance"); err != nil {
		return nil, err
	}
	return &State{X: mat.NewVecDense(len(x), append([]float64(nil), x...)), P: cloneSym(covariance)}, nil
}

// UpdateResult contains the posterior state and update diagnostics.
type UpdateResult struct {
	State                *State
	Gain                 *mat.Dense
	InnovationCovariance *mat.SymDense
	NIS                  float64
}

// Predict applies x(predicted)=F*x and P(predicted)=F*P*F^T+Q.
func Predict(state *State, transition mat.Matrix, processNoise mat.Symmetric) (*State, error) {
	n, err := validateState(state)
	if err != nil {
		return nil, err
	}
	if nilInterface(transition) || nilInterface(processNoise) {
		return nil, fmt.Errorf("sequential: %w: prediction matrices must not be nil", adjust.ErrInvalidProblem)
	}
	fr, fc := transition.Dims()
	if fr != n || fc != n || processNoise.SymmetricDim() != n {
		return nil, fmt.Errorf("sequential: %w: prediction matrix dimensions do not match the state", adjust.ErrInvalidProblem)
	}
	if err := validateFiniteMatrix(transition, "transition matrix"); err != nil {
		return nil, err
	}
	if err := validateCovariance(processNoise, n, "process noise"); err != nil {
		return nil, err
	}
	x := mat.NewVecDense(n, nil)
	x.MulVec(transition, state.X)
	var fp, covariance mat.Dense
	fp.Mul(transition, state.P)
	covariance.Mul(&fp, transition.T())
	for r := 0; r < n; r++ {
		for c := 0; c < n; c++ {
			covariance.Set(r, c, covariance.At(r, c)+processNoise.At(r, c))
		}
	}
	return &State{X: x, P: denseToSym(&covariance)}, nil
}

// Update applies innovation z-h(x), design matrix H, and covariance R using a
// Cholesky solve and Joseph-form covariance update.
func Update(state *State, innovation mat.Vector, design mat.Matrix, measurementCovariance mat.Symmetric) (*UpdateResult, error) {
	n, err := validateState(state)
	if err != nil {
		return nil, err
	}
	if nilInterface(innovation) || nilInterface(design) || nilInterface(measurementCovariance) {
		return nil, fmt.Errorf("sequential: %w: update inputs must not be nil", adjust.ErrInvalidProblem)
	}
	m, columns := design.Dims()
	if m <= 0 || columns != n || innovation.Len() != m || measurementCovariance.SymmetricDim() != m {
		return nil, fmt.Errorf("sequential: %w: update matrix dimensions do not match", adjust.ErrInvalidProblem)
	}
	if err := validateFiniteVector(innovation, "innovation"); err != nil {
		return nil, err
	}
	if err := validateFiniteMatrix(design, "design matrix"); err != nil {
		return nil, err
	}
	if err := validateCovariance(measurementCovariance, m, "measurement covariance"); err != nil {
		return nil, err
	}

	var hp, innovationDense mat.Dense
	hp.Mul(design, state.P)
	innovationDense.Mul(&hp, design.T())
	innovationCovariance := mat.NewSymDense(m, nil)
	for r := 0; r < m; r++ {
		for c := 0; c <= r; c++ {
			value := (innovationDense.At(r, c) + innovationDense.At(c, r)) / 2
			innovationCovariance.SetSym(r, c, value+measurementCovariance.At(r, c))
		}
	}
	var chol mat.Cholesky
	if !chol.Factorize(innovationCovariance) {
		return nil, fmt.Errorf("sequential: %w: innovation covariance", adjust.ErrNotPositiveDefinite)
	}

	var pht, solved mat.Dense
	pht.Mul(state.P, design.T())
	if err := chol.SolveTo(&solved, pht.T()); err != nil {
		return nil, fmt.Errorf("sequential: solve Kalman gain: %w", err)
	}
	gain := mat.DenseCopyOf(solved.T())
	correction := mat.NewVecDense(n, nil)
	correction.MulVec(gain, innovation)
	x := mat.NewVecDense(n, nil)
	x.AddVec(state.X, correction)

	identityMinusKH := mat.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		identityMinusKH.Set(i, i, 1)
	}
	var kh mat.Dense
	kh.Mul(gain, design)
	identityMinusKH.Sub(identityMinusKH, &kh)
	var ap, joseph, kr, noiseTerm mat.Dense
	ap.Mul(identityMinusKH, state.P)
	joseph.Mul(&ap, identityMinusKH.T())
	kr.Mul(gain, measurementCovariance)
	noiseTerm.Mul(&kr, gain.T())
	joseph.Add(&joseph, &noiseTerm)

	weightedInnovation := mat.NewVecDense(m, nil)
	if err := chol.SolveVecTo(weightedInnovation, innovation); err != nil {
		return nil, fmt.Errorf("sequential: solve normalized innovation: %w", err)
	}
	return &UpdateResult{
		State: &State{X: x, P: denseToSym(&joseph)},
		Gain:  gain, InnovationCovariance: innovationCovariance,
		NIS: mat.Dot(innovation, weightedInnovation),
	}, nil
}

func nilInterface(value any) bool {
	if value == nil {
		return true
	}
	reflected := reflect.ValueOf(value)
	switch reflected.Kind() {
	case reflect.Chan, reflect.Func, reflect.Interface, reflect.Map, reflect.Pointer, reflect.Slice:
		return reflected.IsNil()
	default:
		return false
	}
}

func validateState(state *State) (int, error) {
	if state == nil || state.X == nil || state.P == nil || state.X.Len() == 0 || state.P.SymmetricDim() != state.X.Len() {
		return 0, fmt.Errorf("sequential: %w: invalid state", adjust.ErrInvalidProblem)
	}
	if err := validateFiniteVector(state.X, "state"); err != nil {
		return 0, err
	}
	if err := validateCovariance(state.P, state.X.Len(), "state covariance"); err != nil {
		return 0, err
	}
	return state.X.Len(), nil
}

func validateFiniteVector(vector mat.Vector, name string) error {
	for i := 0; i < vector.Len(); i++ {
		value := vector.AtVec(i)
		if math.IsNaN(value) || math.IsInf(value, 0) {
			return fmt.Errorf("sequential: %w: %s contains a non-finite value", adjust.ErrInvalidProblem, name)
		}
	}
	return nil
}

func validateFiniteMatrix(matrix mat.Matrix, name string) error {
	rows, columns := matrix.Dims()
	for row := 0; row < rows; row++ {
		for column := 0; column < columns; column++ {
			value := matrix.At(row, column)
			if math.IsNaN(value) || math.IsInf(value, 0) {
				return fmt.Errorf("sequential: %w: %s contains a non-finite value", adjust.ErrInvalidProblem, name)
			}
		}
	}
	return nil
}

func validateCovariance(covariance mat.Symmetric, dimension int, name string) error {
	if covariance == nil || covariance.SymmetricDim() != dimension {
		return fmt.Errorf("sequential: %w: %s has invalid dimensions", adjust.ErrInvalidProblem, name)
	}
	if err := validateFiniteMatrix(covariance, name); err != nil {
		return err
	}
	var eigen mat.EigenSym
	if !eigen.Factorize(covariance, false) {
		return fmt.Errorf("sequential: %w: cannot decompose %s", adjust.ErrInvalidProblem, name)
	}
	scale := 1.0
	for row := 0; row < dimension; row++ {
		for column := 0; column < dimension; column++ {
			scale = math.Max(scale, math.Abs(covariance.At(row, column)))
		}
	}
	tolerance := 1e-12 * scale
	for _, value := range eigen.Values(nil) {
		if value < -tolerance {
			return fmt.Errorf("sequential: %w: %s must be positive semidefinite", adjust.ErrNotPositiveDefinite, name)
		}
	}
	return nil
}

func denseToSym(dense mat.Matrix) *mat.SymDense {
	n, _ := dense.Dims()
	symmetric := mat.NewSymDense(n, nil)
	for r := 0; r < n; r++ {
		for c := 0; c <= r; c++ {
			value := (dense.At(r, c) + dense.At(c, r)) / 2
			if r == c && value < 0 && value > -1e-14 {
				value = 0
			}
			symmetric.SetSym(r, c, value)
		}
	}
	return symmetric
}

func cloneSym(src *mat.SymDense) *mat.SymDense {
	if src == nil {
		return nil
	}
	return denseToSym(src)
}
