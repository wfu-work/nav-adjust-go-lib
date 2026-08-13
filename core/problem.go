// Package core defines the solver-independent linear adjustment model used by the
// adjustment library.
//
// The package uses the linearized observation convention
//
//	w = H*dx + e
//
// where w is an observed-minus-computed misclosure, H is the design matrix,
// dx is the parameter correction, and e is observation noise.
package core

import "fmt"

// Term is one non-zero coefficient in an observation or constraint row.
type Term struct {
	Parameter   int
	Coefficient float64
}

// T constructs a sparse design-matrix term.
func T(parameter int, coefficient float64) Term {
	return Term{Parameter: parameter, Coefficient: coefficient}
}

// Parameter provides a stable name for one element of the correction vector.
// Name is diagnostic metadata and does not affect the numerical solution.
type Parameter struct {
	Name string
}

// Equation represents one linearized observation equation.
//
// Variance is used when the equation does not belong to a CovarianceBlock. It
// must be positive. Equations in a block take their complete stochastic model
// from that block, so their Variance may be left as zero.
type Equation struct {
	ID         string
	Group      string
	Terms      []Term
	Misclosure float64
	Variance   float64
}

// CovarianceBlock describes correlated observations. RowIndexes refer to
// Problem.Equations and Covariance contains a row-major square matrix.
//
// An equation may belong to at most one covariance block.
type CovarianceBlock struct {
	ID         string
	RowIndexes []int
	Covariance []float64
}

// Constraint represents an exact or stochastic linear parameter constraint.
//
// Exact constraints have Variance == 0. A positive variance creates a soft
// constraint that participates in the stochastic adjustment.
type Constraint struct {
	ID       string
	Terms    []Term
	Value    float64
	Variance float64
}

// ExactConstraint creates C*dx=value as an exact constraint.
func ExactConstraint(id string, value float64, terms ...Term) Constraint {
	return Constraint{ID: id, Terms: terms, Value: value}
}

// SoftConstraint creates C*dx=value with the given variance.
func SoftConstraint(id string, value, variance float64, terms ...Term) Constraint {
	return Constraint{ID: id, Terms: terms, Value: value, Variance: variance}
}

// FixCorrection fixes one element of the correction vector exactly.
func FixCorrection(parameter int, value float64) Constraint {
	return ExactConstraint(fmt.Sprintf("fix:%d", parameter), value, T(parameter, 1))
}

// Problem is a complete linearized least-squares adjustment.
//
// ParameterCount is required. Parameters is optional, but when present it must
// contain ParameterCount entries.
type Problem struct {
	ParameterCount   int
	Parameters       []Parameter
	Equations        []Equation
	CovarianceBlocks []CovarianceBlock
	Constraints      []Constraint
}

// NewProblem creates an empty adjustment problem with optional parameter names.
func NewProblem(parameterNames ...string) Problem {
	problem := Problem{ParameterCount: len(parameterNames)}
	if len(parameterNames) == 0 {
		return problem
	}
	problem.Parameters = make([]Parameter, len(parameterNames))
	for i, name := range parameterNames {
		problem.Parameters[i] = Parameter{Name: name}
	}
	return problem
}

// AddEquation appends an observation equation and returns its row index.
func (p *Problem) AddEquation(e Equation) int {
	p.Equations = append(p.Equations, e)
	return len(p.Equations) - 1
}

// AddConstraint appends a parameter constraint.
func (p *Problem) AddConstraint(c Constraint) {
	p.Constraints = append(p.Constraints, c)
}

// ParameterName returns a stable diagnostic name for a parameter.
func (p Problem) ParameterName(index int) string {
	if index >= 0 && index < len(p.Parameters) && p.Parameters[index].Name != "" {
		return p.Parameters[index].Name
	}
	return fmt.Sprintf("x[%d]", index)
}
