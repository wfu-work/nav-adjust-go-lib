package mock

import (
	"context"
	"fmt"
	"math"
	"runtime"
	"sort"
	"time"

	"github.com/wfu-work/nav-adjust-go-lib/batch"
	"github.com/wfu-work/nav-adjust-go-lib/core"
	"github.com/wfu-work/nav-adjust-go-lib/model"
	"github.com/wfu-work/nav-adjust-go-lib/network"
	"github.com/wfu-work/nav-adjust-go-lib/nonlinear"
)

// Check is one explicit numerical or diagnostic assertion.
type Check struct {
	Name     string
	Expected string
	Actual   string
	Passed   bool
}

// Metric is one display-ready report metric.
type Metric struct {
	Label string
	Value string
	Unit  string
}

// StationComparison compares one estimated station with its known truth.
type StationComparison struct {
	ID       string
	Truth    model.ENU
	Adjusted model.ENU
	Error    float64
}

// ResidualSummary contains one baseline residual ordered by vector norm.
type ResidualSummary struct {
	ID     string
	Group  string
	Norm   float64
	Weight float64
}

// ScenarioResult is a complete report record for one validation case.
type ScenarioResult struct {
	ID          string
	Name        string
	Category    string
	Description string
	Passed      bool
	Duration    time.Duration
	Input       string
	Solver      string
	Metrics     []Metric
	Checks      []Check
	Stations    []StationComparison
	Residuals   []ResidualSummary
	Warnings    []string
	Error       string
}

// CategorySummary aggregates report outcomes by scenario category.
type CategorySummary struct {
	Name   string
	Total  int
	Passed int
}

// SuiteReport is the complete deterministic validation run.
type SuiteReport struct {
	Title       string
	GeneratedAt time.Time
	GoVersion   string
	Duration    time.Duration
	Passed      bool
	Total       int
	PassCount   int
	FailCount   int
	Categories  []CategorySummary
	Scenarios   []ScenarioResult
}

// RunAll executes all network, robust batch, and nonlinear fixtures.
func RunAll(ctx context.Context) SuiteReport {
	started := time.Now()
	report := SuiteReport{
		Title: "NAV Adjust 数值验证报告", GeneratedAt: started, GoVersion: runtime.Version(), Passed: true,
	}
	for _, scenario := range NetworkScenarios() {
		if ctx.Err() != nil {
			report.Scenarios = append(report.Scenarios, canceledScenario(scenario.ID, scenario.Name, scenario.Category, scenario.Description, ctx.Err()))
			break
		}
		report.Scenarios = append(report.Scenarios, runNetworkScenario(ctx, scenario))
	}
	if ctx.Err() == nil {
		report.Scenarios = append(report.Scenarios, runBatchHuberScenario(ctx))
	}
	if ctx.Err() == nil {
		report.Scenarios = append(report.Scenarios, runNonlinearScenario(ctx))
	}
	report.Duration = time.Since(started)
	report.Total = len(report.Scenarios)
	category := make(map[string]*CategorySummary)
	for _, scenario := range report.Scenarios {
		if scenario.Passed {
			report.PassCount++
		} else {
			report.FailCount++
			report.Passed = false
		}
		entry := category[scenario.Category]
		if entry == nil {
			entry = &CategorySummary{Name: scenario.Category}
			category[scenario.Category] = entry
		}
		entry.Total++
		if scenario.Passed {
			entry.Passed++
		}
	}
	for _, entry := range category {
		report.Categories = append(report.Categories, *entry)
	}
	sort.Slice(report.Categories, func(i, j int) bool { return report.Categories[i].Name < report.Categories[j].Name })
	return report
}

func runNetworkScenario(ctx context.Context, scenario NetworkScenario) ScenarioResult {
	started := time.Now()
	record := ScenarioResult{
		ID: scenario.ID, Name: scenario.Name, Category: scenario.Category, Description: scenario.Description,
		Input: fmt.Sprintf("%d 站 / %d 基线 / %d 坐标先验", len(scenario.Problem.Stations), len(scenario.Problem.Baselines), len(scenario.Problem.Priors)),
	}
	if err := network.Validate(scenario.Problem, &scenario.Options); err != nil {
		record.Error = fmt.Sprintf("输入校验失败: %v", err)
		record.Duration = time.Since(started)
		return record
	}
	result, err := network.SolveContext(ctx, scenario.Problem, &scenario.Options)
	record.Duration = time.Since(started)
	if err != nil {
		record.Error = err.Error()
		return record
	}
	record.Solver = result.Diagnostics.Solver
	record.Warnings = append([]string(nil), result.Warnings...)

	comparisons := make([]StationComparison, 0, len(result.Stations))
	squaredError := 0.0
	maxError := 0.0
	for _, station := range result.Stations {
		truth := scenario.Truth[station.ID]
		errorNorm := distanceENU(station.Position, truth)
		squaredError += errorNorm * errorNorm
		maxError = math.Max(maxError, errorNorm)
		comparisons = append(comparisons, StationComparison{ID: station.ID, Truth: truth, Adjusted: station.Position, Error: errorNorm})
	}
	rmse := math.Sqrt(squaredError / math.Max(1, float64(len(comparisons))))
	record.Stations = sampleStations(comparisons, 12)
	record.Residuals = largestResiduals(result.Baselines, 8)
	record.Metrics = []Metric{
		{Label: "坐标 RMSE", Value: formatNumber(rmse, 6), Unit: "m"},
		{Label: "最大三维误差", Value: formatNumber(maxError, 6), Unit: "m"},
		{Label: "目标函数", Value: formatNumber(result.Diagnostics.Objective, 6)},
		{Label: "单位权中误差", Value: formatNumber(result.Diagnostics.Sigma0, 6)},
		{Label: "自由度", Value: fmt.Sprint(result.Diagnostics.DegreesOfFreedom)},
		{Label: "求解迭代", Value: fmt.Sprint(result.Diagnostics.SolverIterations)},
	}
	addCheck(&record, "坐标真值误差", fmt.Sprintf("<= %.6g m", scenario.Expectation.MaxPositionError), fmt.Sprintf("%.6g m", maxError), maxError <= scenario.Expectation.MaxPositionError)
	addCheck(&record, "参数秩", "rank = parameter count", fmt.Sprintf("%d / %d", result.Diagnostics.Rank, result.Diagnostics.ParameterCount), result.Diagnostics.Rank == result.Diagnostics.ParameterCount)
	addCheck(&record, "迭代状态", "converged", fmt.Sprint(result.Diagnostics.Converged), result.Diagnostics.Converged)
	addCheck(&record, "求解器路径", scenario.Expectation.Solver, result.Diagnostics.Solver, result.Diagnostics.Solver == scenario.Expectation.Solver)
	if scenario.Expectation.Preconditioner != "" {
		addCheck(&record, "预条件器", scenario.Expectation.Preconditioner, result.Diagnostics.SolverPreconditioner, result.Diagnostics.SolverPreconditioner == scenario.Expectation.Preconditioner)
	}
	addCheck(&record, "基准模式", string(scenario.Expectation.Datum), string(result.Diagnostics.DatumMode), result.Diagnostics.DatumMode == scenario.Expectation.Datum)
	checkCovarianceMode(&record, result.Diagnostics, scenario.Expectation.CovarianceMode)
	if scenario.Expectation.DownweightedBaselineID != "" {
		baseline, found := findBaseline(result.Baselines, scenario.Expectation.DownweightedBaselineID)
		actual := "not found"
		passed := false
		if found {
			actual = fmt.Sprintf("weight=%.6g, downweighted=%v", baseline.Weight, baseline.Downweighted)
			passed = baseline.Downweighted && baseline.Weight < 0.5
		}
		addCheck(&record, "粗差基线降权", "downweighted with weight < 0.5", actual, passed)
	}
	if scenario.Expectation.LargerVarianceGroup != "" {
		large, largeFound := varianceScale(result.VarianceComponents, scenario.Expectation.LargerVarianceGroup)
		small, smallFound := varianceScale(result.VarianceComponents, scenario.Expectation.SmallerVarianceGroup)
		ratio := large / small
		actual := fmt.Sprintf("%s=%.6g, %s=%.6g, ratio=%.3g", scenario.Expectation.LargerVarianceGroup, large, scenario.Expectation.SmallerVarianceGroup, small, ratio)
		addCheck(&record, "方差分量区分", fmt.Sprintf("ratio >= %.3g", scenario.Expectation.MinimumVarianceRatio), actual, largeFound && smallFound && ratio >= scenario.Expectation.MinimumVarianceRatio)
		addCheck(&record, "方差分量迭代", "converged", fmt.Sprint(result.Diagnostics.VarianceComponentsConverged), result.Diagnostics.VarianceComponentsConverged)
	}
	record.Passed = checksPassed(record.Checks)
	return record
}

func runBatchHuberScenario(ctx context.Context) ScenarioResult {
	started := time.Now()
	record := ScenarioResult{
		ID: "batch-huber", Name: "通用 Huber 线性平差", Category: "稳健估计",
		Description: "同一标量参数包含五个正常观测和一个大粗差，验证 batch.SolveHuber 的逐观测重加权。",
		Input:       "1 参数 / 6 独立观测",
	}
	problem := core.Problem{ParameterCount: 1}
	for i, value := range []float64{10.02, 9.99, 10.01, 10.03, 9.98, 18} {
		problem.Equations = append(problem.Equations, core.Equation{
			ID: fmt.Sprintf("obs-%d", i+1), Terms: []core.Term{core.T(0, 1)}, Misclosure: value, Variance: 0.01 * 0.01,
		})
	}
	result, err := batch.SolveHuberContext(ctx, problem, &batch.HuberOptions{K: 1.5, MaxIterations: 40, Tolerance: 1e-6, MinWeight: 0.001})
	record.Duration = time.Since(started)
	if err != nil {
		record.Error = err.Error()
		return record
	}
	estimate := result.Adjustment.Delta[0]
	errorValue := math.Abs(estimate - 10)
	record.Solver = result.Adjustment.Method
	record.Metrics = []Metric{
		{Label: "估计值", Value: formatNumber(estimate, 6)},
		{Label: "绝对误差", Value: formatNumber(errorValue, 6)},
		{Label: "粗差权", Value: formatNumber(result.Weights[5], 6)},
		{Label: "迭代次数", Value: fmt.Sprint(result.Iterations)},
	}
	addCheck(&record, "参数真值误差", "<= 0.03", formatNumber(errorValue, 6), errorValue <= 0.03)
	addCheck(&record, "粗差观测降权", "weight < 0.1", formatNumber(result.Weights[5], 6), result.Weights[5] < 0.1)
	addCheck(&record, "Huber 迭代", "converged", fmt.Sprint(result.Converged), result.Converged)
	record.Passed = checksPassed(record.Checks)
	return record
}

func runNonlinearScenario(ctx context.Context) ScenarioResult {
	started := time.Now()
	record := ScenarioResult{
		ID: "nonlinear-spp", Name: "非线性伪距定位", Category: "非线性平差",
		Description: "六颗卫星的合成伪距通过 Gauss-Newton 反复线性化，联合估计 ECEF 坐标和接收机钟差。",
		Input:       "4 参数 / 6 伪距观测",
	}
	truth := [4]float64{-2267800, 5009300, 3221000, 75000}
	positions := [][3]float64{
		{15600e3, 7540e3, 20140e3}, {18760e3, 2750e3, 18610e3},
		{17610e3, 14630e3, 13480e3}, {19170e3, 610e3, 18390e3},
		{17800e3, -15200e3, 15100e3}, {-13400e3, 21300e3, 9000e3},
	}
	noise := []float64{1.2, -0.8, 0.4, -1, 0.7, -0.2}
	pseudorange := make([]float64, len(positions))
	for i, position := range positions {
		pseudorange[i] = distance3(truth[:3], position[:]) + truth[3] + noise[i]
	}
	linearModel := nonlinear.ModelFunc(func(state []float64) (core.Problem, error) {
		problem := core.NewProblem("ecef-x", "ecef-y", "ecef-z", "clock-metres")
		for i, position := range positions {
			dx, dy, dz := state[0]-position[0], state[1]-position[1], state[2]-position[2]
			rho := math.Sqrt(dx*dx + dy*dy + dz*dz)
			problem.AddEquation(core.Equation{
				ID: fmt.Sprintf("G%02d-C1C", i+1), Group: "GPS-C1C",
				Terms:      []core.Term{core.T(0, dx/rho), core.T(1, dy/rho), core.T(2, dz/rho), core.T(3, 1)},
				Misclosure: pseudorange[i] - (rho + state[3]), Variance: 4,
			})
		}
		return problem, nil
	})
	result, err := nonlinear.SolveContext(ctx, []float64{-2266800, 5008300, 3221500, 0}, linearModel, nil)
	record.Duration = time.Since(started)
	if err != nil {
		record.Error = err.Error()
		return record
	}
	positionError := distance3(result.State[:3], truth[:3])
	clockError := math.Abs(result.State[3] - truth[3])
	record.Solver = result.Adjustment.Method
	record.Metrics = []Metric{
		{Label: "ECEF 三维误差", Value: formatNumber(positionError, 6), Unit: "m"},
		{Label: "钟差误差", Value: formatNumber(clockError, 6), Unit: "m"},
		{Label: "Gauss-Newton 迭代", Value: fmt.Sprint(result.Iterations)},
		{Label: "单位权中误差", Value: formatNumber(result.Adjustment.Sigma0, 6)},
	}
	addCheck(&record, "ECEF 真值误差", "<= 1 m", fmt.Sprintf("%.6g m", positionError), positionError <= 1)
	addCheck(&record, "钟差真值误差", "<= 1 m", fmt.Sprintf("%.6g m", clockError), clockError <= 1)
	addCheck(&record, "非线性收敛", "converged within 6 iterations", fmt.Sprintf("converged=%v, iterations=%d", result.Converged, result.Iterations), result.Converged && result.Iterations <= 6)
	record.Passed = checksPassed(record.Checks)
	return record
}

func checkCovarianceMode(record *ScenarioResult, diagnostics model.NetworkDiagnostics, expected model.CovarianceMode) {
	actual := string(diagnostics.CovarianceMode)
	passed := diagnostics.CovarianceMode == expected
	switch expected {
	case model.CovarianceFull:
		passed = passed && diagnostics.FullCovarianceAvailable && diagnostics.StationCovarianceAvailable && diagnostics.ResidualDiagnosticsAvailable
	case model.CovarianceStationBlocks:
		passed = passed && !diagnostics.FullCovarianceAvailable && diagnostics.StationCovarianceAvailable && diagnostics.ResidualDiagnosticsAvailable
	case model.CovarianceNone:
		passed = passed && !diagnostics.FullCovarianceAvailable && !diagnostics.StationCovarianceAvailable && !diagnostics.ResidualDiagnosticsAvailable
	}
	addCheck(record, "协方差输出模式", string(expected), actual, passed)
}

func addCheck(record *ScenarioResult, name, expected, actual string, passed bool) {
	record.Checks = append(record.Checks, Check{Name: name, Expected: expected, Actual: actual, Passed: passed})
}

func checksPassed(checks []Check) bool {
	if len(checks) == 0 {
		return false
	}
	for _, check := range checks {
		if !check.Passed {
			return false
		}
	}
	return true
}

func findBaseline(baselines []model.BaselineResult, id string) (model.BaselineResult, bool) {
	for _, baseline := range baselines {
		if baseline.ID == id {
			return baseline, true
		}
	}
	return model.BaselineResult{}, false
}

func varianceScale(components []model.VarianceComponentResult, group string) (float64, bool) {
	for _, component := range components {
		if component.Group == group {
			return component.Scale, true
		}
	}
	return 0, false
}

func largestResiduals(baselines []model.BaselineResult, limit int) []ResidualSummary {
	result := make([]ResidualSummary, len(baselines))
	for i, baseline := range baselines {
		result[i] = ResidualSummary{ID: baseline.ID, Group: baseline.Group, Norm: normENU(baseline.Residual), Weight: baseline.Weight}
	}
	sort.SliceStable(result, func(i, j int) bool { return result[i].Norm > result[j].Norm })
	if len(result) > limit {
		result = result[:limit]
	}
	return result
}

func sampleStations(stations []StationComparison, limit int) []StationComparison {
	if len(stations) <= limit {
		return stations
	}
	front := limit - 3
	result := append([]StationComparison(nil), stations[:front]...)
	return append(result, stations[len(stations)-3:]...)
}

func canceledScenario(id, name, category, description string, err error) ScenarioResult {
	return ScenarioResult{ID: id, Name: name, Category: category, Description: description, Error: err.Error()}
}

func distanceENU(left, right model.ENU) float64 {
	return normENU(subtractENU(left, right))
}

func normENU(value model.ENU) float64 {
	return math.Sqrt(value.East*value.East + value.North*value.North + value.Up*value.Up)
}

func distance3(left, right []float64) float64 {
	dx, dy, dz := left[0]-right[0], left[1]-right[1], left[2]-right[2]
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

func formatNumber(value float64, precision int) string {
	return fmt.Sprintf("%.*f", precision, value)
}
