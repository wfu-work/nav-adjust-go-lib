package network

import "github.com/wfu-work/nav-adjust-go-lib/batch"

type networkOptions struct {
	confidence     float64
	conditionLimit float64
	batch          batch.Options
	robust         *robustNetworkOptions
	variance       *varianceComponentNetworkOptions
	covariance     CovarianceMode
	datum          DatumMode
}

type robustNetworkOptions struct {
	threshold     float64
	maxIterations int
	tolerance     float64
	minWeight     float64
}

type varianceComponentNetworkOptions struct {
	maxIterations     int
	tolerance         float64
	minScale          float64
	maxScale          float64
	minimumRedundancy float64
}

func normalizeNetworkOptions(options *ENUNetworkOptions) (networkOptions, error) {
	public := ENUNetworkOptions{}
	if options != nil {
		public = *options
	}
	result := networkOptions{}
	result.datum = public.Datum
	if result.datum == "" {
		result.datum = DatumExternal
	}
	if result.datum != DatumExternal && result.datum != DatumFreeCentroid {
		return networkOptions{}, validationError(ErrUnsupportedMethod, "options.datum", string(result.datum), "must be external or free-centroid")
	}
	result.covariance = public.Covariance
	if result.covariance == "" {
		result.covariance = CovarianceFull
	}
	if result.covariance != CovarianceFull && result.covariance != CovarianceStationBlocks && result.covariance != CovarianceNone {
		return networkOptions{}, validationError(nil, "options.covariance", string(result.covariance), "must be full, station-blocks, or none")
	}
	result.confidence = public.Confidence
	if result.confidence == 0 {
		result.confidence = 0.95
	}
	if result.confidence <= 0 || result.confidence >= 1 || !finite(result.confidence) {
		return networkOptions{}, validationError(nil, "options.confidence", "", "must be between zero and one")
	}
	result.conditionLimit = public.ConditionWarnLimit
	if result.conditionLimit == 0 {
		result.conditionLimit = 1e12
	}
	if result.conditionLimit <= 0 || !finite(result.conditionLimit) {
		return networkOptions{}, validationError(nil, "options.condition_warn_limit", "", "must be positive")
	}
	result.batch = batch.Options{
		SymmetryTolerance:   public.SymmetryTolerance,
		MinVariance:         public.MinimumVariance,
		DenseThreshold:      public.Solver.DenseThreshold,
		MaxIterations:       public.Solver.MaxIterations,
		RelativeTolerance:   public.Solver.RelativeTolerance,
		AbsoluteTolerance:   public.Solver.AbsoluteTolerance,
		PreconditionerShift: public.Solver.PreconditionerShift,
	}
	switch public.Solver.Method {
	case "", SolverDense:
		result.batch.Solver = batch.SolverDense
	case SolverAuto:
		result.batch.Solver = batch.SolverAuto
	case SolverPCG:
		result.batch.Solver = batch.SolverPCG
	default:
		return networkOptions{}, validationError(ErrUnsupportedMethod, "options.solver.method", string(public.Solver.Method), "only dense, auto, and pcg are supported")
	}
	switch public.Solver.Preconditioner {
	case "":
		if public.Solver.Method == SolverAuto {
			result.batch.Preconditioner = batch.PreconditionerIC0
		} else {
			result.batch.Preconditioner = batch.PreconditionerBlockJacobi
			result.batch.PreconditionerBlockSize = 3
		}
	case PreconditionerBlockJacobi:
		result.batch.Preconditioner = batch.PreconditionerBlockJacobi
		result.batch.PreconditionerBlockSize = 3
	case PreconditionerJacobi:
		result.batch.Preconditioner = batch.PreconditionerJacobi
	case PreconditionerIC0:
		result.batch.Preconditioner = batch.PreconditionerIC0
	default:
		return networkOptions{}, validationError(ErrUnsupportedMethod, "options.solver.preconditioner", string(public.Solver.Preconditioner), "only jacobi, block-jacobi, and ic0 are supported")
	}
	if result.covariance == CovarianceFull {
		result.batch.Covariance = batch.CovarianceFull
	} else {
		result.batch.Covariance = batch.CovarianceNone
	}
	if public.SymmetryTolerance < 0 || !finite(public.SymmetryTolerance) {
		return networkOptions{}, validationError(nil, "options.symmetry_tolerance", "", "must be non-negative")
	}
	if public.MinimumVariance < 0 || !finite(public.MinimumVariance) {
		return networkOptions{}, validationError(nil, "options.minimum_variance", "", "must be non-negative")
	}
	if public.Solver.MaxIterations < 0 {
		return networkOptions{}, validationError(nil, "options.solver.max_iterations", "", "must be non-negative")
	}
	if public.Solver.DenseThreshold < 0 {
		return networkOptions{}, validationError(nil, "options.solver.dense_threshold", "", "must be non-negative")
	}
	if public.Solver.PreconditionerShift < 0 || !finite(public.Solver.PreconditionerShift) {
		return networkOptions{}, validationError(nil, "options.solver.preconditioner_shift", "", "must be finite and non-negative")
	}
	if public.Solver.RelativeTolerance < 0 || !finite(public.Solver.RelativeTolerance) {
		return networkOptions{}, validationError(nil, "options.solver.relative_tolerance", "", "must be non-negative")
	}
	if public.Solver.AbsoluteTolerance < 0 || !finite(public.Solver.AbsoluteTolerance) {
		return networkOptions{}, validationError(nil, "options.solver.absolute_tolerance", "", "must be non-negative")
	}
	if public.VarianceComponents != nil {
		variance := &varianceComponentNetworkOptions{
			maxIterations:     public.VarianceComponents.MaxIterations,
			tolerance:         public.VarianceComponents.Tolerance,
			minScale:          public.VarianceComponents.MinScale,
			maxScale:          public.VarianceComponents.MaxScale,
			minimumRedundancy: public.VarianceComponents.MinimumRedundancy,
		}
		if variance.maxIterations == 0 {
			variance.maxIterations = 10
		}
		if variance.tolerance == 0 {
			variance.tolerance = 1e-3
		}
		if variance.minScale == 0 {
			variance.minScale = 1e-4
		}
		if variance.maxScale == 0 {
			variance.maxScale = 1e4
		}
		if variance.minimumRedundancy == 0 {
			variance.minimumRedundancy = 1e-6
		}
		if variance.maxIterations <= 0 {
			return networkOptions{}, validationError(nil, "options.variance_components.max_iterations", "", "must be positive")
		}
		if variance.tolerance <= 0 || !finite(variance.tolerance) {
			return networkOptions{}, validationError(nil, "options.variance_components.tolerance", "", "must be positive")
		}
		if variance.minScale <= 0 || !finite(variance.minScale) {
			return networkOptions{}, validationError(nil, "options.variance_components.min_scale", "", "must be positive")
		}
		if variance.maxScale < variance.minScale || !finite(variance.maxScale) {
			return networkOptions{}, validationError(nil, "options.variance_components.max_scale", "", "must be no smaller than min_scale")
		}
		if variance.minimumRedundancy <= 0 || !finite(variance.minimumRedundancy) {
			return networkOptions{}, validationError(nil, "options.variance_components.minimum_redundancy", "", "must be positive")
		}
		result.variance = variance
	}
	if public.Robust == nil {
		return result, nil
	}
	method := public.Robust.Method
	if method == "" {
		method = RobustHuber
	}
	if method != RobustHuber {
		return networkOptions{}, validationError(ErrUnsupportedMethod, "options.robust.method", string(method), "only huber is supported")
	}
	robust := &robustNetworkOptions{
		threshold:     public.Robust.Threshold,
		maxIterations: public.Robust.MaxIterations,
		tolerance:     public.Robust.Tolerance,
		minWeight:     public.Robust.MinWeight,
	}
	if robust.threshold == 0 {
		robust.threshold = 2.5
	}
	if robust.maxIterations == 0 {
		robust.maxIterations = 10
	}
	if robust.tolerance == 0 {
		robust.tolerance = 1e-3
	}
	if robust.minWeight == 0 {
		robust.minWeight = 0.05
	}
	if robust.threshold <= 0 || !finite(robust.threshold) {
		return networkOptions{}, validationError(nil, "options.robust.threshold", "", "must be positive")
	}
	if robust.maxIterations <= 0 {
		return networkOptions{}, validationError(nil, "options.robust.max_iterations", "", "must be positive")
	}
	if robust.tolerance <= 0 || !finite(robust.tolerance) {
		return networkOptions{}, validationError(nil, "options.robust.tolerance", "", "must be positive")
	}
	if robust.minWeight <= 0 || robust.minWeight > 1 || !finite(robust.minWeight) {
		return networkOptions{}, validationError(nil, "options.robust.min_weight", "", "must be in (0, 1]")
	}
	result.robust = robust
	return result, nil
}
