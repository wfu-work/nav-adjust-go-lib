package core

import "errors"

var (
	// ErrInvalidProblem indicates inconsistent dimensions, indices, or numbers.
	ErrInvalidProblem = errors.New("adjust: invalid problem")
	// ErrRankDeficient indicates that observations and constraints do not define
	// all requested parameters.
	ErrRankDeficient = errors.New("adjust: rank deficient")
	// ErrNotPositiveDefinite indicates an invalid covariance or normal matrix.
	ErrNotPositiveDefinite = errors.New("adjust: matrix is not positive definite")
	// ErrNotConverged indicates that an iterative numerical solver exhausted its
	// iteration limit before reaching the requested residual tolerance.
	ErrNotConverged = errors.New("adjust: iterative solver did not converge")
)
