package batch

import (
	"context"
	"errors"
	"fmt"
	"math"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
	"github.com/wfu-work/nav-adjust-go-lib/internal/sparse"
	"gonum.org/v1/gonum/mat"
)

// Solution contains the ordinary adjustment result and retains the numerical
// factorization needed for explicit, on-demand covariance blocks. It is useful
// to higher-level packages that need local covariance without an n-by-n matrix.
type Solution struct {
	Result   *adjust.Result
	provider covarianceProvider
	size     int
}

// FormalCovarianceBlock returns the unscaled inverse-normal submatrix in the
// order supplied by indexes. It performs iterative column solves in PCG mode
// and does not retain complete inverse columns.
func (solution *Solution) FormalCovarianceBlock(indexes []int) (*mat.SymDense, error) {
	return solution.FormalCovarianceBlockContext(context.Background(), indexes)
}

// FormalCovarianceBlockContext is FormalCovarianceBlock with cancellation for
// iterative covariance-column solves.
func (solution *Solution) FormalCovarianceBlockContext(ctx context.Context, indexes []int) (*mat.SymDense, error) {
	if solution == nil || solution.Result == nil || solution.provider == nil {
		return nil, fmt.Errorf("batch: covariance is unavailable")
	}
	for _, index := range indexes {
		if index < 0 || index >= solution.size {
			return nil, invalid("covariance parameter index %d is out of range", index)
		}
	}
	return covarianceBlock(ctx, solution.provider, solution.size, indexes)
}

// FormalCovarianceValues returns selected inverse-normal entries. Each pair is
// [row, column], and the output follows the input order. Columns are solved at
// most once per call, so callers can request a sparse set of station and edge
// covariance entries without materializing the complete inverse.
func (solution *Solution) FormalCovarianceValues(pairs [][2]int) ([]float64, error) {
	return solution.FormalCovarianceValuesContext(context.Background(), pairs)
}

// FormalCovarianceValuesContext is FormalCovarianceValues with cancellation
// for iterative covariance-column solves.
func (solution *Solution) FormalCovarianceValuesContext(ctx context.Context, pairs [][2]int) ([]float64, error) {
	if solution == nil || solution.Result == nil || solution.provider == nil {
		return nil, fmt.Errorf("batch: covariance is unavailable")
	}
	byColumn := make(map[int][]int)
	for i, pair := range pairs {
		if pair[0] < 0 || pair[0] >= solution.size || pair[1] < 0 || pair[1] >= solution.size {
			return nil, invalid("covariance parameter pair [%d,%d] is out of range", pair[0], pair[1])
		}
		byColumn[pair[1]] = append(byColumn[pair[1]], i)
	}
	values := make([]float64, len(pairs))
	for column, positions := range byColumn {
		if err := contextError(ctx); err != nil {
			return nil, err
		}
		data, err := solution.provider.column(ctx, column)
		if err != nil {
			return nil, err
		}
		for _, position := range positions {
			values[position] = data[pairs[position][0]]
		}
	}
	return values, nil
}

// Solve performs a linear weighted least-squares adjustment. Zero-value
// options retain the historical dense solve with complete covariance.
func Solve(problem adjust.Problem, options *Options) (*adjust.Result, error) {
	return SolveContext(context.Background(), problem, options)
}

// SolveContext is Solve with cooperative cancellation.
func SolveContext(ctx context.Context, problem adjust.Problem, options *Options) (*adjust.Result, error) {
	solution, err := SolveDetailedContext(ctx, problem, options)
	if err != nil {
		return nil, err
	}
	return solution.Result, nil
}

// SolveDetailed solves a problem and also permits on-demand covariance block
// queries. CovarianceNone omits complete inverse construction from Result.
func SolveDetailed(problem adjust.Problem, options *Options) (*Solution, error) {
	return SolveDetailedContext(context.Background(), problem, options)
}

// SolveDetailedContext is SolveDetailed with cooperative cancellation.
func SolveDetailedContext(ctx context.Context, problem adjust.Problem, options *Options) (*Solution, error) {
	if err := contextError(ctx); err != nil {
		return nil, err
	}
	opts, err := normalizeOptions(options)
	if err != nil {
		return nil, err
	}
	if err := contextError(ctx); err != nil {
		return nil, err
	}
	compiled, err := compileProblem(problem, opts)
	if err != nil {
		return nil, err
	}
	solver := resolveSolver(opts, compiled.parameters)

	var delta []float64
	var provider covarianceProvider
	method := ""
	condition := 0.0
	conditionAvailable := false
	iterations := 0
	relativeResidual := 0.0

	if solver == SolverPCG {
		normal, buildErr := buildSparseNormal(ctx, compiled)
		if buildErr != nil {
			return nil, buildErr
		}
		pcgOptions := sparse.PCGOptions{
			MaxIterations: opts.MaxIterations, RelativeTolerance: opts.RelativeTolerance,
			AbsoluteTolerance: opts.AbsoluteTolerance,
		}
		if compiled.exactCount == 0 {
			pcgOptions.Preconditioner, err = makePreconditioner(normal.n, opts, nil)
			if err != nil {
				return nil, mapPCGError("preconditioner", err)
			}
			delta, iterations, relativeResidual, err = solveSparse(ctx, normal, pcgOptions)
			if err != nil {
				return nil, err
			}
			provider = &pcgCovarianceProvider{matrix: normal.n, options: pcgOptions}
			method = "sparse-pcg"
		} else {
			delta, provider, iterations, relativeResidual, err = solveProjectedSparse(ctx, normal, compiled, opts, pcgOptions)
			if err != nil {
				return nil, err
			}
			method = "sparse-projected-pcg"
		}
	} else {
		normal, buildErr := buildNormal(ctx, compiled)
		if buildErr != nil {
			return nil, buildErr
		}
		var denseDelta *mat.VecDense
		denseDelta, provider, method, condition, err = solveDense(ctx, normal, compiled)
		if err != nil {
			return nil, err
		}
		delta = append([]float64(nil), denseDelta.RawVector().Data...)
		conditionAvailable = true
	}

	var formal *mat.SymDense
	if opts.Covariance == CovarianceFull {
		formal, err = fullCovariance(ctx, provider, compiled.parameters)
		if err != nil {
			return nil, err
		}
	}
	result, err := makeResult(compiled, delta, formal, method, condition)
	if err != nil {
		return nil, err
	}
	result.ConditionAvailable = conditionAvailable
	result.SolverIterations = iterations
	result.SolverRelativeResidual = relativeResidual
	if pcgOptions, ok := providerPCGOptions(provider); ok && pcgOptions.Preconditioner != nil {
		result.SolverPreconditioner = pcgOptions.Preconditioner.Name()
	}
	return &Solution{Result: result, provider: provider, size: compiled.parameters}, nil
}

type covarianceProvider interface {
	column(ctx context.Context, index int) ([]float64, error)
}

type choleskyCovarianceProvider struct {
	factor *mat.Cholesky
	size   int
}

func (provider *choleskyCovarianceProvider) column(ctx context.Context, index int) ([]float64, error) {
	if err := contextError(ctx); err != nil {
		return nil, err
	}
	rhs := mat.NewVecDense(provider.size, nil)
	rhs.SetVec(index, 1)
	result := mat.NewVecDense(provider.size, nil)
	if err := provider.factor.SolveVecTo(result, rhs); err != nil {
		return nil, fmt.Errorf("batch: covariance column %d: %w", index, err)
	}
	return append([]float64(nil), result.RawVector().Data...), nil
}

type kktCovarianceProvider struct {
	factor     *mat.LU
	parameters int
	size       int
}

func (provider *kktCovarianceProvider) column(ctx context.Context, index int) ([]float64, error) {
	if err := contextError(ctx); err != nil {
		return nil, err
	}
	rhs := mat.NewDense(provider.size, 1, nil)
	rhs.Set(index, 0, 1)
	var result mat.Dense
	if err := provider.factor.SolveTo(&result, false, rhs); err != nil {
		return nil, fmt.Errorf("batch: constrained covariance column %d: %w", index, err)
	}
	column := make([]float64, provider.parameters)
	for i := range column {
		column[i] = result.At(i, 0)
	}
	return column, nil
}

type pcgCovarianceProvider struct {
	matrix  *sparse.Matrix
	options sparse.PCGOptions
}

func (provider *pcgCovarianceProvider) column(ctx context.Context, index int) ([]float64, error) {
	rhs := make([]float64, provider.matrix.Size())
	rhs[index] = 1
	column, _, err := sparse.SolvePCGContext(ctx, provider.matrix, rhs, provider.options)
	if err != nil {
		return nil, mapPCGError("covariance query", err)
	}
	return column, nil
}

func solveSparse(ctx context.Context, normal sparseNormalEquation, options sparse.PCGOptions) ([]float64, int, float64, error) {
	delta, info, err := sparse.SolvePCGContext(ctx, normal.n, normal.u, options)
	if err != nil {
		return nil, info.Iterations, info.RelativeResidual, mapPCGError("normal equation", err)
	}
	return delta, info.Iterations, info.RelativeResidual, nil
}

func mapPCGError(context string, err error) error {
	switch {
	case errors.Is(err, sparse.ErrNotConverged):
		return fmt.Errorf("%w: %s: %v", adjust.ErrNotConverged, context, err)
	case errors.Is(err, sparse.ErrBreakdown):
		return fmt.Errorf("%w: %s: %v", adjust.ErrRankDeficient, context, err)
	default:
		return fmt.Errorf("batch: sparse %s: %w", context, err)
	}
}

func solveDense(ctx context.Context, normal normalEquation, compiled compiledProblem) (*mat.VecDense, covarianceProvider, string, float64, error) {
	if err := contextError(ctx); err != nil {
		return nil, nil, "", 0, err
	}
	nx := compiled.parameters
	if compiled.exactCount == 0 {
		factor := &mat.Cholesky{}
		if !factor.Factorize(normal.n) {
			return nil, nil, "cholesky", math.Inf(1), fmt.Errorf("%w: normal equation", adjust.ErrRankDeficient)
		}
		if err := contextError(ctx); err != nil {
			return nil, nil, "cholesky", factor.Cond(), err
		}
		delta := mat.NewVecDense(nx, nil)
		if err := factor.SolveVecTo(delta, normal.u); err != nil {
			return nil, nil, "cholesky", math.Inf(1), fmt.Errorf("batch: solve normal equation: %w", err)
		}
		if err := contextError(ctx); err != nil {
			return nil, nil, "cholesky", factor.Cond(), err
		}
		return delta, &choleskyCovarianceProvider{factor: factor, size: nx}, "cholesky", factor.Cond(), nil
	}

	size := nx + compiled.exactCount
	kkt := mat.NewDense(size, size, nil)
	for r := 0; r < nx; r++ {
		for c := 0; c < nx; c++ {
			kkt.Set(r, c, normal.n.At(r, c))
		}
	}
	for r := 0; r < compiled.exactCount; r++ {
		for _, term := range compiled.exact[r].terms {
			kkt.Set(term.Parameter, nx+r, term.Coefficient)
			kkt.Set(nx+r, term.Parameter, term.Coefficient)
		}
	}
	rhs := mat.NewDense(size, 1, nil)
	for i := 0; i < nx; i++ {
		rhs.Set(i, 0, normal.u.AtVec(i))
	}
	for i := 0; i < compiled.exactCount; i++ {
		rhs.Set(nx+i, 0, compiled.exact[i].value)
	}
	factor := &mat.LU{}
	factor.Factorize(kkt)
	condition := factor.Cond()
	if err := contextError(ctx); err != nil {
		return nil, nil, "kkt", condition, err
	}
	if math.IsInf(condition, 1) || condition > 1e16 {
		return nil, nil, "kkt", condition, fmt.Errorf("%w: constrained normal equation has condition number %.6g", adjust.ErrRankDeficient, condition)
	}
	var rawSolution mat.Dense
	if err := factor.SolveTo(&rawSolution, false, rhs); err != nil {
		if errors.Is(err, mat.ErrSingular) {
			return nil, nil, "kkt", condition, fmt.Errorf("%w: constrained normal equation: %v", adjust.ErrRankDeficient, err)
		}
		return nil, nil, "kkt", condition, fmt.Errorf("batch: solve constrained normal equation: %w", err)
	}
	if err := contextError(ctx); err != nil {
		return nil, nil, "kkt", condition, err
	}
	delta := mat.NewVecDense(nx, nil)
	for i := range nx {
		delta.SetVec(i, rawSolution.At(i, 0))
	}
	return delta, &kktCovarianceProvider{factor: factor, parameters: nx, size: size}, "kkt", condition, nil
}

func fullCovariance(ctx context.Context, provider covarianceProvider, size int) (*mat.SymDense, error) {
	if err := contextError(ctx); err != nil {
		return nil, err
	}
	switch provider := provider.(type) {
	case *choleskyCovarianceProvider:
		result := mat.NewSymDense(size, nil)
		provider.factor.InverseTo(result)
		return result, nil
	case *kktCovarianceProvider:
		identity := mat.NewDense(provider.size, provider.size, nil)
		for i := range provider.size {
			identity.Set(i, i, 1)
		}
		var inverse mat.Dense
		if err := provider.factor.SolveTo(&inverse, false, identity); err != nil {
			return nil, fmt.Errorf("batch: covariance of constrained adjustment: %w", err)
		}
		result := mat.NewSymDense(size, nil)
		for row := range size {
			for column := 0; column <= row; column++ {
				result.SetSym(row, column, (inverse.At(row, column)+inverse.At(column, row))/2)
			}
		}
		return result, nil
	}
	indexes := make([]int, size)
	for i := range indexes {
		indexes[i] = i
	}
	return covarianceBlock(ctx, provider, size, indexes)
}

func covarianceBlock(ctx context.Context, provider covarianceProvider, parameterCount int, indexes []int) (*mat.SymDense, error) {
	columns := make([][]float64, len(indexes))
	for i, index := range indexes {
		if err := contextError(ctx); err != nil {
			return nil, err
		}
		column, err := provider.column(ctx, index)
		if err != nil {
			return nil, err
		}
		if len(column) != parameterCount {
			return nil, fmt.Errorf("batch: covariance provider returned an invalid column")
		}
		columns[i] = column
	}
	result := mat.NewSymDense(len(indexes), nil)
	for r := range indexes {
		for c := 0; c <= r; c++ {
			value := (columns[c][indexes[r]] + columns[r][indexes[c]]) / 2
			result.SetSym(r, c, value)
		}
	}
	return result, nil
}

func contextError(ctx context.Context) error {
	if ctx == nil {
		return fmt.Errorf("batch: nil context")
	}
	return ctx.Err()
}

func makePreconditioner(matrix *sparse.Matrix, options Options, project sparse.ProjectFunc) (sparse.Preconditioner, error) {
	switch options.Preconditioner {
	case PreconditionerBlockJacobi:
		return sparse.NewBlockJacobiPreconditioner(matrix, options.PreconditionerBlockSize, project)
	case PreconditionerIC0:
		return sparse.NewIncompleteCholeskyPreconditioner(matrix, options.PreconditionerShift)
	default:
		return sparse.NewJacobiPreconditioner(matrix, project)
	}
}

func resolveSolver(options Options, parameterCount int) SolverMethod {
	if options.Solver != SolverAuto {
		return options.Solver
	}
	if options.Covariance == CovarianceFull || parameterCount <= options.DenseThreshold {
		return SolverDense
	}
	return SolverPCG
}

func providerPCGOptions(provider covarianceProvider) (sparse.PCGOptions, bool) {
	switch provider := provider.(type) {
	case *pcgCovarianceProvider:
		return provider.options, true
	case *projectedPCGCovarianceProvider:
		return provider.options, true
	default:
		return sparse.PCGOptions{}, false
	}
}

func makeResult(compiled compiledProblem, delta []float64, formal *mat.SymDense, method string, condition float64) (*adjust.Result, error) {
	result := &adjust.Result{
		Delta: append([]float64(nil), delta...), FormalCovariance: formal,
		Rank: compiled.parameters, Condition: condition, Method: method,
		ExactConstraintCount: compiled.exactCount,
		CovarianceAvailable:  formal != nil, ResidualDiagnosticsAvailable: formal != nil,
	}
	for _, row := range compiled.rows {
		if row.soft {
			result.SoftConstraintCount++
		} else {
			result.ObservationCount++
		}
	}
	result.DOF = len(compiled.rows) + compiled.exactCount - compiled.parameters

	values := make([]float64, len(compiled.rows))
	for i, row := range compiled.rows {
		values[i] = row.w - dotTerms(row.terms, result.Delta)
	}
	objective, err := weightedObjective(compiled, values)
	if err != nil {
		return nil, err
	}
	result.Objective = objective
	result.Sigma0 = 1
	if result.DOF > 0 {
		result.Sigma0 = math.Sqrt(math.Max(0, objective/float64(result.DOF)))
	}
	if formal != nil {
		result.Covariance = scaleSym(formal, result.Sigma0*result.Sigma0)
	}

	var residualVariances, redundancies []float64
	if formal != nil {
		residualVariances, redundancies, err = residualDiagnostics(compiled, formal)
		if err != nil {
			return nil, err
		}
	}
	for i, row := range compiled.rows {
		qvv, redundancy, standardized := 0.0, 0.0, 0.0
		if formal != nil {
			qvv, redundancy = residualVariances[i], redundancies[i]
			if qvv > 0 && result.Sigma0 > 0 {
				standardized = values[i] / (result.Sigma0 * math.Sqrt(qvv))
			}
		}
		if row.soft {
			result.ConstraintResiduals = append(result.ConstraintResiduals, adjust.ConstraintResidual{
				ID: row.id, Value: values[i], Variance: qvv, Standardized: standardized,
			})
			continue
		}
		result.Residuals = append(result.Residuals, adjust.Residual{
			ID: row.id, Group: row.group, Value: values[i], Variance: qvv,
			Redundancy: redundancy, Standardized: standardized,
		})
	}
	return result, nil
}

func residualDiagnostics(compiled compiledProblem, covariance mat.Symmetric) ([]float64, []float64, error) {
	variances := make([]float64, len(compiled.rows))
	redundancies := make([]float64, len(compiled.rows))
	for blockIndex := range compiled.blocks {
		block := &compiled.blocks[blockIndex]
		n := len(block.rowIndexes)
		residualCovariance := mat.NewSymDense(n, nil)
		for r, rowR := range block.rowIndexes {
			for c := 0; c <= r; c++ {
				value := block.covariance.At(r, c)
				rowC := block.rowIndexes[c]
				residualCovariance.SetSym(r, c, value-sparseCrossQuadratic(compiled.rows[rowR].terms, covariance, compiled.rows[rowC].terms))
			}
		}
		var weightedResidualCovariance mat.Dense
		if err := block.cholesky.SolveTo(&weightedResidualCovariance, residualCovariance); err != nil {
			return nil, nil, fmt.Errorf("batch: residual diagnostics for block %q: %w", block.id, err)
		}
		for local, row := range block.rowIndexes {
			variances[row] = nonnegativeRoundoff(residualCovariance.At(local, local), compiled.rows[row].variance)
			redundancies[row] = clamp(weightedResidualCovariance.At(local, local), 0, 1)
		}
	}
	for i, row := range compiled.rows {
		if compiled.blocked[i] {
			continue
		}
		variances[i] = nonnegativeRoundoff(row.variance-sparseCrossQuadratic(row.terms, covariance, row.terms), row.variance)
		redundancies[i] = clamp(variances[i]/row.variance, 0, 1)
	}
	return variances, redundancies, nil
}

func nonnegativeRoundoff(value, scale float64) float64 {
	if value < 0 && value > -1e-10*math.Max(1, scale) {
		return 0
	}
	return value
}

func weightedObjective(compiled compiledProblem, residuals []float64) (float64, error) {
	objective := 0.0
	for blockIndex := range compiled.blocks {
		block := &compiled.blocks[blockIndex]
		n := len(block.rowIndexes)
		v := mat.NewVecDense(n, nil)
		for r, row := range block.rowIndexes {
			v.SetVec(r, residuals[row])
		}
		weighted := mat.NewVecDense(n, nil)
		if err := block.cholesky.SolveVecTo(weighted, v); err != nil {
			return 0, err
		}
		objective += mat.Dot(v, weighted)
	}
	for i, row := range compiled.rows {
		if !compiled.blocked[i] {
			objective += residuals[i] * residuals[i] / row.variance
		}
	}
	return objective, nil
}

func sparseCrossQuadratic(left []adjust.Term, covariance mat.Symmetric, right []adjust.Term) float64 {
	value := 0.0
	for _, leftTerm := range left {
		for _, rightTerm := range right {
			value += leftTerm.Coefficient * covariance.At(leftTerm.Parameter, rightTerm.Parameter) * rightTerm.Coefficient
		}
	}
	return value
}

func dotTerms(terms []adjust.Term, values []float64) float64 {
	value := 0.0
	for _, term := range terms {
		value += term.Coefficient * values[term.Parameter]
	}
	return value
}

func scaleSym(src *mat.SymDense, scale float64) *mat.SymDense {
	n := src.SymmetricDim()
	dst := mat.NewSymDense(n, nil)
	for r := 0; r < n; r++ {
		for c := 0; c <= r; c++ {
			dst.SetSym(r, c, src.At(r, c)*scale)
		}
	}
	return dst
}

func clamp(value, low, high float64) float64 {
	return math.Max(low, math.Min(high, value))
}
