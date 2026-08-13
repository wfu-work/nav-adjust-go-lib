package batch

import (
	"fmt"
	"math"
	"sort"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
	"gonum.org/v1/gonum/mat"
)

type stochasticRow struct {
	id       string
	group    string
	terms    []adjust.Term
	w        float64
	variance float64
	soft     bool
}

type stochasticBlock struct {
	id         string
	rowIndexes []int
	parameters []int
	covariance *mat.SymDense
	cholesky   mat.Cholesky
}

type exactRow struct {
	terms []adjust.Term
	value float64
}

type compiledProblem struct {
	rows       []stochasticRow
	blocks     []stochasticBlock
	blocked    []bool
	exact      []exactRow
	exactCount int
	parameters int
}

func compileProblem(problem adjust.Problem, options Options) (compiledProblem, error) {
	if problem.ParameterCount <= 0 {
		return compiledProblem{}, invalid("parameter count must be positive")
	}
	if len(problem.Parameters) != 0 && len(problem.Parameters) != problem.ParameterCount {
		return compiledProblem{}, invalid("parameter metadata count %d does not match parameter count %d", len(problem.Parameters), problem.ParameterCount)
	}
	if len(problem.Equations) == 0 && len(problem.Constraints) == 0 {
		return compiledProblem{}, invalid("problem has no equations or constraints")
	}

	compiled := compiledProblem{parameters: problem.ParameterCount}
	compiled.rows = make([]stochasticRow, 0, len(problem.Equations)+len(problem.Constraints))
	for i, equation := range problem.Equations {
		terms, err := canonicalTerms(problem.ParameterCount, equation.Terms)
		if err != nil {
			return compiledProblem{}, invalid("equation %q: %v", rowID(equation.ID, i), err)
		}
		if !finite(equation.Misclosure) {
			return compiledProblem{}, invalid("equation %q has a non-finite misclosure", rowID(equation.ID, i))
		}
		compiled.rows = append(compiled.rows, stochasticRow{
			id: equation.ID, group: equation.Group, terms: terms,
			w: equation.Misclosure, variance: equation.Variance,
		})
	}

	blocked := make([]bool, len(problem.Equations))
	compiled.blocks = make([]stochasticBlock, len(problem.CovarianceBlocks))
	for blockIndex, block := range problem.CovarianceBlocks {
		n := len(block.RowIndexes)
		if n == 0 || len(block.Covariance) != n*n {
			return compiledProblem{}, invalid("covariance block %q must contain an n-by-n row-major matrix", block.ID)
		}
		for local, row := range block.RowIndexes {
			if row < 0 || row >= len(problem.Equations) {
				return compiledProblem{}, invalid("covariance block %q row index %d is out of range", block.ID, row)
			}
			if blocked[row] {
				return compiledProblem{}, invalid("equation row %d belongs to more than one covariance block", row)
			}
			blocked[row] = true
			diagonal := block.Covariance[local*n+local]
			if !finite(diagonal) || diagonal < options.MinVariance {
				return compiledProblem{}, invalid("covariance block %q has invalid variance at row %d", block.ID, local)
			}
			compiled.rows[row].variance = diagonal
		}
		for r := 0; r < n; r++ {
			for c := 0; c < n; c++ {
				value := block.Covariance[r*n+c]
				if !finite(value) {
					return compiledProblem{}, invalid("covariance block %q contains a non-finite value", block.ID)
				}
				if math.Abs(value-block.Covariance[c*n+r]) > options.SymmetryTolerance {
					return compiledProblem{}, invalid("covariance block %q is not symmetric", block.ID)
				}
			}
		}
		covariance := mat.NewSymDense(n, nil)
		for r := 0; r < n; r++ {
			for c := 0; c <= r; c++ {
				covariance.SetSym(r, c, block.Covariance[r*n+c])
			}
		}
		var chol mat.Cholesky
		if !chol.Factorize(covariance) {
			return compiledProblem{}, fmt.Errorf("%w: covariance block %q", adjust.ErrNotPositiveDefinite, block.ID)
		}
		condition := chol.Cond()
		if math.IsInf(condition, 1) || condition > 1e16 {
			return compiledProblem{}, fmt.Errorf(
				"%w: covariance block %q has condition number %.6g",
				adjust.ErrNotPositiveDefinite, block.ID, condition,
			)
		}
		compiled.blocks[blockIndex] = stochasticBlock{
			id:         block.ID,
			rowIndexes: append([]int(nil), block.RowIndexes...),
			parameters: blockParameters(compiled.rows, block.RowIndexes),
			covariance: covariance,
			cholesky:   chol,
		}
	}
	for row, equation := range problem.Equations {
		if blocked[row] {
			continue
		}
		if !finite(equation.Variance) || equation.Variance < options.MinVariance {
			return compiledProblem{}, invalid("equation %q must have a positive variance", rowID(equation.ID, row))
		}
	}

	exact := make([]exactRow, 0, len(problem.Constraints))
	for i, constraint := range problem.Constraints {
		terms, err := canonicalTerms(problem.ParameterCount, constraint.Terms)
		if err != nil {
			return compiledProblem{}, invalid("constraint %q: %v", rowID(constraint.ID, i), err)
		}
		if !finite(constraint.Value) || !finite(constraint.Variance) || constraint.Variance < 0 {
			return compiledProblem{}, invalid("constraint %q has an invalid value or variance", rowID(constraint.ID, i))
		}
		if constraint.Variance == 0 {
			exact = append(exact, exactRow{terms: terms, value: constraint.Value})
			continue
		}
		if constraint.Variance < options.MinVariance {
			return compiledProblem{}, invalid("soft constraint %q variance is too small; use an exact constraint", rowID(constraint.ID, i))
		}
		compiled.rows = append(compiled.rows, stochasticRow{
			id: constraint.ID, group: "constraint", terms: terms,
			w: constraint.Value, variance: constraint.Variance, soft: true,
		})
	}
	compiled.exactCount = len(exact)
	compiled.exact = exact
	compiled.blocked = make([]bool, len(compiled.rows))
	copy(compiled.blocked, blocked)
	return compiled, nil
}

func canonicalTerms(parameterCount int, terms []adjust.Term) ([]adjust.Term, error) {
	if len(terms) == 0 {
		return nil, fmt.Errorf("row has no terms")
	}
	coefficients := make(map[int]float64, len(terms))
	for _, term := range terms {
		if term.Parameter < 0 || term.Parameter >= parameterCount {
			return nil, fmt.Errorf("parameter index %d is out of range", term.Parameter)
		}
		if !finite(term.Coefficient) {
			return nil, fmt.Errorf("coefficient for parameter %d is not finite", term.Parameter)
		}
		coefficient := coefficients[term.Parameter] + term.Coefficient
		if !finite(coefficient) {
			return nil, fmt.Errorf("combined coefficient for parameter %d is not finite", term.Parameter)
		}
		coefficients[term.Parameter] = coefficient
	}
	parameters := make([]int, 0, len(coefficients))
	for parameter, coefficient := range coefficients {
		if coefficient != 0 {
			parameters = append(parameters, parameter)
		}
	}
	if len(parameters) == 0 {
		return nil, fmt.Errorf("row has no non-zero coefficients")
	}
	sort.Ints(parameters)
	canonical := make([]adjust.Term, len(parameters))
	for i, parameter := range parameters {
		canonical[i] = adjust.T(parameter, coefficients[parameter])
	}
	return canonical, nil
}

func blockParameters(rows []stochasticRow, rowIndexes []int) []int {
	set := make(map[int]struct{})
	for _, rowIndex := range rowIndexes {
		for _, term := range rows[rowIndex].terms {
			set[term.Parameter] = struct{}{}
		}
	}
	parameters := make([]int, 0, len(set))
	for parameter := range set {
		parameters = append(parameters, parameter)
	}
	sort.Ints(parameters)
	return parameters
}

func finite(value float64) bool {
	return !math.IsNaN(value) && !math.IsInf(value, 0)
}

func rowID(id string, index int) string {
	if id != "" {
		return id
	}
	return fmt.Sprintf("row:%d", index)
}

func invalid(format string, args ...any) error {
	return fmt.Errorf("%w: %s", adjust.ErrInvalidProblem, fmt.Sprintf(format, args...))
}
