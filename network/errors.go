package network

import (
	"errors"
	"fmt"

	"github.com/wfu-work/nav-adjust-go-lib/core"
)

var (
	// ErrInvalidProblem indicates invalid input fields or topology.
	ErrInvalidProblem = core.ErrInvalidProblem
	// ErrRankDeficient indicates that the network does not determine every free station.
	ErrRankDeficient = core.ErrRankDeficient
	// ErrNotPositiveDefinite indicates an invalid covariance or normal matrix.
	ErrNotPositiveDefinite = core.ErrNotPositiveDefinite
	// ErrNotConverged indicates that the sparse iterative solver reached its
	// configured iteration limit.
	ErrNotConverged = core.ErrNotConverged

	// ErrInvalidCovariance indicates an invalid observation covariance matrix.
	ErrInvalidCovariance = errors.New("adjust: invalid covariance")
	// ErrDuplicateStation indicates repeated station IDs.
	ErrDuplicateStation = errors.New("adjust: duplicate station")
	// ErrDuplicateBaseline indicates repeated baseline IDs.
	ErrDuplicateBaseline = errors.New("adjust: duplicate baseline")
	// ErrDuplicatePrior indicates repeated position-prior IDs.
	ErrDuplicatePrior = errors.New("adjust: duplicate position prior")
	// ErrUnknownStation indicates that a referenced station does not exist.
	ErrUnknownStation = errors.New("adjust: unknown station")
	// ErrDisconnectedNetwork indicates a component without a usable datum or
	// relative observation.
	ErrDisconnectedNetwork = errors.New("adjust: disconnected network component has no fixed station or position prior")
	// ErrUnsupportedMethod indicates an unknown robust or solver method.
	ErrUnsupportedMethod = errors.New("adjust: unsupported method")
	// ErrInsufficientRedundancy indicates that at least one baseline group does
	// not have enough independent residual information to estimate its scale.
	ErrInsufficientRedundancy = core.ErrInsufficientRedundancy
)

// ValidationError identifies a public input field and record that failed
// validation. errors.Is also matches ErrInvalidProblem and Kind.
type ValidationError struct {
	Kind    error
	Field   string
	ID      string
	Message string
}

func (err *ValidationError) Error() string {
	location := err.Field
	if err.ID != "" {
		location += "[" + err.ID + "]"
	}
	if err.Message == "" {
		return fmt.Sprintf("adjust: invalid %s", location)
	}
	return fmt.Sprintf("adjust: invalid %s: %s", location, err.Message)
}

// Is supports both broad and specific errors.Is checks.
func (err *ValidationError) Is(target error) bool {
	return target == ErrInvalidProblem || errors.Is(err.Kind, target)
}

// Unwrap returns the specific validation cause, when present.
func (err *ValidationError) Unwrap() error {
	return err.Kind
}

func validationError(kind error, field, id, message string) error {
	return &ValidationError{Kind: kind, Field: field, ID: id, Message: message}
}
