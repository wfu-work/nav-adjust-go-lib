package model

// AdjustedStation is one fixed or adjusted station in the common ENU frame.
type AdjustedStation struct {
	ID               string            `json:"id"`
	Name             string            `json:"name,omitempty"`
	Fixed            bool              `json:"fixed"`
	Position         ENU               `json:"position"`
	StdDev           ENU               `json:"stddev"`
	FormalCovariance Matrix3           `json:"formal_covariance"`
	Covariance       Matrix3           `json:"covariance"`
	Metadata         map[string]string `json:"metadata,omitempty"`
}

// BaselineResult contains the observed and adjusted vector and diagnostics.
// Residual follows Observed - Adjusted. Weight is one block-wise robust weight
// for the complete ENU vector and equals one without robust adjustment.
type BaselineResult struct {
	ID                       string            `json:"id"`
	From                     string            `json:"from"`
	To                       string            `json:"to"`
	Group                    string            `json:"group,omitempty"`
	Observed                 ENU               `json:"observed"`
	Adjusted                 ENU               `json:"adjusted"`
	Residual                 ENU               `json:"residual"`
	ResidualStdDev           ENU               `json:"residual_stddev"`
	Standardized             ENU               `json:"standardized"`
	Redundancy               ENU               `json:"redundancy"`
	FormalResidualCovariance Matrix3           `json:"formal_residual_covariance"`
	ResidualCovariance       Matrix3           `json:"residual_covariance"`
	Weight                   float64           `json:"weight"`
	Downweighted             bool              `json:"downweighted"`
	VarianceScale            float64           `json:"variance_scale"`
	Metadata                 map[string]string `json:"metadata,omitempty"`
}

// PositionPriorResult contains the adjusted station coordinate and diagnostics
// for one stochastic control coordinate. Residual follows Observed - Adjusted.
type PositionPriorResult struct {
	ID                       string            `json:"id"`
	StationID                string            `json:"station_id"`
	Observed                 ENU               `json:"observed"`
	Adjusted                 ENU               `json:"adjusted"`
	Residual                 ENU               `json:"residual"`
	ResidualStdDev           ENU               `json:"residual_stddev"`
	Standardized             ENU               `json:"standardized"`
	Redundancy               ENU               `json:"redundancy"`
	FormalResidualCovariance Matrix3           `json:"formal_residual_covariance"`
	ResidualCovariance       Matrix3           `json:"residual_covariance"`
	Metadata                 map[string]string `json:"metadata,omitempty"`
}

// GlobalTestResult is a two-sided chi-square model test.
type GlobalTestResult struct {
	Statistic  float64 `json:"statistic"`
	DOF        int     `json:"degrees_of_freedom"`
	Confidence float64 `json:"confidence"`
	Lower      float64 `json:"lower"`
	Upper      float64 `json:"upper"`
	Passed     bool    `json:"passed"`
}

// VarianceComponentResult describes the estimated covariance multiplier for
// one baseline group. StdDevScale is sqrt(Scale).
type VarianceComponentResult struct {
	Group            string  `json:"group"`
	Scale            float64 `json:"scale"`
	StdDevScale      float64 `json:"stddev_scale"`
	BaselineCount    int     `json:"baseline_count"`
	ObservationCount int     `json:"observation_count"`
	Objective        float64 `json:"objective"`
	Redundancy       float64 `json:"redundancy"`
}

// NetworkDiagnostics summarizes the adjustment and numerical solution.
type NetworkDiagnostics struct {
	StationCount                 int               `json:"station_count"`
	FixedStationCount            int               `json:"fixed_station_count"`
	FreeStationCount             int               `json:"free_station_count"`
	BaselineCount                int               `json:"baseline_count"`
	PriorCount                   int               `json:"prior_count"`
	ObservationCount             int               `json:"observation_count"`
	ParameterCount               int               `json:"parameter_count"`
	Rank                         int               `json:"rank"`
	DegreesOfFreedom             int               `json:"degrees_of_freedom"`
	Objective                    float64           `json:"objective"`
	Sigma0                       float64           `json:"sigma0"`
	ConditionNumber              float64           `json:"condition_number"`
	ConditionNumberAvailable     bool              `json:"condition_number_available"`
	Solver                       string            `json:"solver"`
	SolverPreconditioner         string            `json:"solver_preconditioner,omitempty"`
	SolverIterations             int               `json:"solver_iterations"`
	SolverRelativeResidual       float64           `json:"solver_relative_residual"`
	CovarianceMode               CovarianceMode    `json:"covariance_mode"`
	FullCovarianceAvailable      bool              `json:"full_covariance_available"`
	StationCovarianceAvailable   bool              `json:"station_covariance_available"`
	ResidualDiagnosticsAvailable bool              `json:"residual_diagnostics_available"`
	DatumMode                    DatumMode         `json:"datum_mode"`
	FreeDatumComponentCount      int               `json:"free_datum_component_count"`
	InternalDatumConstraintCount int               `json:"internal_datum_constraint_count"`
	Iterations                   int               `json:"iterations"`
	Converged                    bool              `json:"converged"`
	VarianceComponentIterations  int               `json:"variance_component_iterations"`
	VarianceComponentsAvailable  bool              `json:"variance_components_available"`
	VarianceComponentsConverged  bool              `json:"variance_components_converged"`
	GlobalTest                   *GlobalTestResult `json:"global_test,omitempty"`
}

// ENUNetworkResult is the complete public output. In CovarianceFull mode both
// covariance matrices follow ParameterKeys and include cross-station terms.
// Reduced modes leave the top-level matrices empty; availability is explicit
// in Diagnostics.
type ENUNetworkResult struct {
	Name               string                    `json:"name,omitempty"`
	Stations           []AdjustedStation         `json:"stations"`
	Baselines          []BaselineResult          `json:"baselines"`
	Priors             []PositionPriorResult     `json:"priors,omitempty"`
	ParameterKeys      []string                  `json:"parameter_keys"`
	FormalCovariance   Matrix                    `json:"formal_covariance"`
	Covariance         Matrix                    `json:"covariance"`
	Diagnostics        NetworkDiagnostics        `json:"diagnostics"`
	VarianceComponents []VarianceComponentResult `json:"variance_components,omitempty"`
	Warnings           []string                  `json:"warnings,omitempty"`
}

// Result is the short name used with network.Solve.
type Result = ENUNetworkResult
