// Package mock provides deterministic adjustment fixtures, truth data, and
// validation helpers for examples, regression tests, and report generation.
package mock

import (
	"fmt"
	"math"

	"github.com/wfu-work/nav-adjust-go-lib/model"
)

// NetworkExpectation describes numerical and diagnostic checks for a fixture.
type NetworkExpectation struct {
	MaxPositionError       float64
	Solver                 string
	Preconditioner         string
	Datum                  model.DatumMode
	DownweightedBaselineID string
	LargerVarianceGroup    string
	SmallerVarianceGroup   string
	MinimumVarianceRatio   float64
	CovarianceMode         model.CovarianceMode
}

// NetworkScenario is one deterministic ENU network and its known truth.
type NetworkScenario struct {
	ID          string
	Name        string
	Category    string
	Description string
	Problem     model.Problem
	Options     model.Options
	Truth       map[string]model.ENU
	Expectation NetworkExpectation
}

// NetworkScenarios returns independent copies of all built-in network cases.
func NetworkScenarios() []NetworkScenario {
	return []NetworkScenario{
		fixedDatumScenario(),
		correlatedCovarianceScenario(),
		stochasticControlScenario(),
		freeNetworkScenario(),
		robustOutlierScenario(),
		varianceComponentScenario(),
		sparseAutoScenario(150),
	}
}

func fixedDatumScenario() NetworkScenario {
	truth := compactTruth()
	problem := networkProblem("fixed-datum", truth, "A")
	addBaselines(&problem, truth, []baselineFixture{
		{"AB", "A", "B", model.ENU{East: 0.004, North: -0.003, Up: 0.002}, ""},
		{"BC", "B", "C", model.ENU{East: -0.002, North: 0.004, Up: -0.001}, ""},
		{"CD", "C", "D", model.ENU{East: 0.003, North: 0.001, Up: -0.002}, ""},
		{"DA", "D", "A", model.ENU{East: -0.004, North: -0.002, Up: 0.001}, ""},
		{"AC", "A", "C", model.ENU{East: 0.001, North: -0.001, Up: 0.001}, ""},
		{"BD", "B", "D", model.ENU{East: -0.001, North: 0.002, Up: -0.001}, ""},
	}, model.Matrix3FromStdDev(0.02, 0.02, 0.03))
	return NetworkScenario{
		ID: "fixed-datum", Name: "固定站闭合网", Category: "经典网平差",
		Description: "一个固定控制站和六条冗余 ENU 基线，验证稠密 Cholesky、闭合约束与完整协方差输出。",
		Problem:     problem, Options: model.Options{Covariance: model.CovarianceFull}, Truth: truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.015, Solver: "cholesky", Datum: model.DatumExternal, CovarianceMode: model.CovarianceFull},
	}
}

func correlatedCovarianceScenario() NetworkScenario {
	truth := compactTruth()
	problem := networkProblem("correlated-covariance", truth, "A")
	covariance := model.Matrix3{Data: [9]float64{
		0.0004, 0.00012, -0.00006,
		0.00012, 0.0005, 0.00008,
		-0.00006, 0.00008, 0.0009,
	}}
	addBaselines(&problem, truth, []baselineFixture{
		{"AB-cor", "A", "B", model.ENU{East: 0.008, North: 0.003, Up: -0.004}, ""},
		{"BC-cor", "B", "C", model.ENU{East: -0.006, North: 0.009, Up: 0.003}, ""},
		{"CD-cor", "C", "D", model.ENU{East: 0.004, North: -0.005, Up: 0.006}, ""},
		{"DA-cor", "D", "A", model.ENU{East: -0.005, North: -0.004, Up: -0.002}, ""},
		{"AC-cor", "A", "C", model.ENU{East: 0.003, North: 0.002, Up: -0.003}, ""},
		{"BD-cor", "B", "D", model.ENU{East: -0.002, North: 0.004, Up: 0.002}, ""},
	}, covariance)
	return NetworkScenario{
		ID: "correlated-covariance", Name: "相关协方差网", Category: "随机模型",
		Description: "每条观测携带含 E/N/U 相关项的完整 3x3 协方差，验证非对角随机模型和精度传播。",
		Problem:     problem, Options: model.Options{Covariance: model.CovarianceFull}, Truth: truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.025, Solver: "cholesky", Datum: model.DatumExternal, CovarianceMode: model.CovarianceFull},
	}
}

func stochasticControlScenario() NetworkScenario {
	truth := compactTruth()
	problem := networkProblem("stochastic-control", truth, "")
	addBaselines(&problem, truth, []baselineFixture{
		{"AB-soft", "A", "B", model.ENU{East: 0.002, North: -0.001}, ""},
		{"BC-soft", "B", "C", model.ENU{North: 0.002, Up: -0.001}, ""},
		{"CD-soft", "C", "D", model.ENU{East: -0.001, Up: 0.001}, ""},
		{"DA-soft", "D", "A", model.ENU{East: -0.001, North: -0.001}, ""},
		{"AC-soft", "A", "C", model.ENU{East: 0.001, North: 0.001}, ""},
	}, model.Matrix3FromStdDev(0.015, 0.015, 0.025))
	problem.Priors = []model.PositionPrior{
		{ID: "prior-A", StationID: "A", Position: addENU(truth["A"], model.ENU{East: 0.01, North: -0.015, Up: 0.005}), Covariance: model.Matrix3FromStdDev(0.05, 0.05, 0.08)},
		{ID: "prior-C", StationID: "C", Position: addENU(truth["C"], model.ENU{East: -0.02, North: 0.01, Up: -0.01}), Covariance: model.Matrix3FromStdDev(0.08, 0.08, 0.1)},
	}
	return NetworkScenario{
		ID: "stochastic-control", Name: "软控制点网", Category: "基准模型",
		Description: "不设置固定站，使用两个带协方差的坐标先验建立外部基准，验证软控制点残差和站点协方差块。",
		Problem:     problem, Options: model.Options{Covariance: model.CovarianceStationBlocks}, Truth: truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.025, Solver: "cholesky", Datum: model.DatumExternal, CovarianceMode: model.CovarianceStationBlocks},
	}
}

func freeNetworkScenario() NetworkScenario {
	truth := map[string]model.ENU{
		"A": {East: -10, North: -4, Up: -1},
		"B": {East: 0, North: -6, Up: 0.5},
		"C": {East: 12, North: -1, Up: 1},
		"D": {East: 8, North: 7, Up: -0.5},
		"E": {East: -10, North: 4, Up: 0},
	}
	problem := networkProblem("free-centroid", truth, "")
	addBaselines(&problem, truth, []baselineFixture{
		{"AB-free", "A", "B", model.ENU{East: 0.001, North: -0.001}, ""},
		{"BC-free", "B", "C", model.ENU{North: 0.001, Up: -0.001}, ""},
		{"CD-free", "C", "D", model.ENU{East: -0.001, Up: 0.001}, ""},
		{"DE-free", "D", "E", model.ENU{East: 0.001, North: 0.001}, ""},
		{"EA-free", "E", "A", model.ENU{East: -0.001, North: -0.001}, ""},
		{"AC-free", "A", "C", model.ENU{East: 0.001, Up: 0.001}, ""},
		{"BD-free", "B", "D", model.ENU{North: -0.001, Up: -0.001}, ""},
	}, model.Matrix3FromStdDev(0.01, 0.01, 0.015))
	return NetworkScenario{
		ID: "free-centroid", Name: "零质心自由网", Category: "基准模型",
		Description: "全网无固定站和先验，以三个精确零质心约束定义内部基准，并强制走稀疏投影 PCG。",
		Problem:     problem,
		Options:     model.Options{Datum: model.DatumFreeCentroid, Covariance: model.CovarianceNone, Solver: model.SolverOptions{Method: model.SolverAuto, DenseThreshold: 1}},
		Truth:       truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.004, Solver: "sparse-projected-pcg", Preconditioner: "ic0", Datum: model.DatumFreeCentroid, CovarianceMode: model.CovarianceNone},
	}
}

func robustOutlierScenario() NetworkScenario {
	truth := compactTruth()
	problem := networkProblem("robust-outlier", truth, "A")
	addBaselines(&problem, truth, []baselineFixture{
		{"AB-good", "A", "B", model.ENU{East: 0.004, North: -0.002}, "survey"},
		{"BC-good", "B", "C", model.ENU{East: -0.003, North: 0.003}, "survey"},
		{"CD-good", "C", "D", model.ENU{East: 0.002, North: -0.002}, "survey"},
		{"DA-good", "D", "A", model.ENU{East: -0.003, North: 0.001}, "survey"},
		{"BD-good", "B", "D", model.ENU{East: 0.002, North: 0.002}, "survey"},
		{"CA-good", "C", "A", model.ENU{East: -0.001, North: -0.002}, "survey"},
		{"gross-AC", "A", "C", model.ENU{East: 3, North: -2, Up: 1}, "survey"},
	}, model.Matrix3FromStdDev(0.02, 0.02, 0.03))
	return NetworkScenario{
		ID: "robust-outlier", Name: "粗差鲁棒网", Category: "稳健估计",
		Description: "在冗余闭合网中注入一条米级粗差，验证整条 ENU 基线的 Mahalanobis Huber 降权。",
		Problem:     problem,
		Options:     model.Options{Covariance: model.CovarianceFull, Robust: &model.RobustOptions{Method: model.RobustHuber, Threshold: 2.5, MaxIterations: 30, Tolerance: 1e-5, MinWeight: 0.001}},
		Truth:       truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.04, Solver: "cholesky", Datum: model.DatumExternal, DownweightedBaselineID: "gross-AC", CovarianceMode: model.CovarianceFull},
	}
}

func varianceComponentScenario() NetworkScenario {
	truth := compactTruth()
	problem := networkProblem("variance-components", truth, "A")
	precise := []baselineFixture{
		{"p-AB", "A", "B", model.ENU{East: 0.003, North: -0.002}, "precise"},
		{"p-BC", "B", "C", model.ENU{East: -0.004, North: 0.003}, "precise"},
		{"p-CD", "C", "D", model.ENU{East: 0.002, North: -0.003}, "precise"},
		{"p-DA", "D", "A", model.ENU{East: -0.002, North: 0.001}, "precise"},
		{"p-AC", "A", "C", model.ENU{East: 0.003, Up: -0.002}, "precise"},
		{"p-BD", "B", "D", model.ENU{North: 0.002, Up: 0.002}, "precise"},
	}
	noisy := []baselineFixture{
		{"n-AB", "A", "B", model.ENU{East: 0.11, North: -0.08, Up: 0.04}, "noisy"},
		{"n-BC", "B", "C", model.ENU{East: -0.09, North: 0.12, Up: -0.05}, "noisy"},
		{"n-CD", "C", "D", model.ENU{East: 0.08, North: -0.1, Up: 0.06}, "noisy"},
		{"n-DA", "D", "A", model.ENU{East: -0.12, North: 0.07, Up: -0.03}, "noisy"},
		{"n-AC", "A", "C", model.ENU{East: 0.1, North: 0.09, Up: -0.07}, "noisy"},
		{"n-BD", "B", "D", model.ENU{East: -0.08, North: 0.11, Up: 0.05}, "noisy"},
	}
	addBaselines(&problem, truth, append(precise, noisy...), model.Matrix3FromStdDev(0.02, 0.02, 0.03))
	return NetworkScenario{
		ID: "variance-components", Name: "分组方差分量", Category: "随机模型",
		Description: "同一网形包含精密组和噪声被低估组，验证组协方差尺度的迭代估计与排序。",
		Problem:     problem,
		Options:     model.Options{Covariance: model.CovarianceStationBlocks, VarianceComponents: &model.VarianceComponentOptions{MaxIterations: 30, Tolerance: 1e-3, MinScale: 1e-4, MaxScale: 1e4, MinimumRedundancy: 1e-6}},
		Truth:       truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.04, Solver: "cholesky", Datum: model.DatumExternal, LargerVarianceGroup: "noisy", SmallerVarianceGroup: "precise", MinimumVarianceRatio: 10, CovarianceMode: model.CovarianceStationBlocks},
	}
}

func sparseAutoScenario(stationCount int) NetworkScenario {
	truth := make(map[string]model.ENU, stationCount)
	problem := model.Problem{Name: "sparse-auto"}
	for i := 0; i < stationCount; i++ {
		id := fmt.Sprintf("S%03d", i)
		position := model.ENU{East: float64(i) * 1.2, North: math.Sin(float64(i)*0.08) * 3, Up: float64(i) * 0.01}
		truth[id] = position
		station := model.Station{ID: id}
		if i == 0 {
			known := position
			station.Fixed = true
			station.KnownENU = &known
		}
		problem.Stations = append(problem.Stations, station)
	}
	fixtures := make([]baselineFixture, 0, stationCount*2)
	for i := 0; i < stationCount-1; i++ {
		fixtures = append(fixtures, baselineFixture{
			ID: fmt.Sprintf("chain-%03d", i), From: fmt.Sprintf("S%03d", i), To: fmt.Sprintf("S%03d", i+1),
			Noise: patternedNoise(i, 0.001), Group: "chain",
		})
	}
	for i := 0; i < stationCount-3; i++ {
		fixtures = append(fixtures, baselineFixture{
			ID: fmt.Sprintf("brace-%03d", i), From: fmt.Sprintf("S%03d", i), To: fmt.Sprintf("S%03d", i+3),
			Noise: patternedNoise(i+2, 0.0015), Group: "brace",
		})
	}
	addBaselines(&problem, truth, fixtures, model.Matrix3FromStdDev(0.02, 0.02, 0.03))
	return NetworkScenario{
		ID: "sparse-auto-150", Name: "150 站稀疏网", Category: "规模与求解器",
		Description: "150 个站点、链边与跨三站加固边，验证 Auto 阈值、CSR 法方程和 IC(0)-PCG。",
		Problem:     problem,
		Options:     model.Options{Covariance: model.CovarianceNone, Solver: model.SolverOptions{Method: model.SolverAuto, DenseThreshold: 300}},
		Truth:       truth,
		Expectation: NetworkExpectation{MaxPositionError: 0.05, Solver: "sparse-pcg", Preconditioner: "ic0", Datum: model.DatumExternal, CovarianceMode: model.CovarianceNone},
	}
}

type baselineFixture struct {
	ID    string
	From  string
	To    string
	Noise model.ENU
	Group string
}

func compactTruth() map[string]model.ENU {
	return map[string]model.ENU{
		"A": {East: 100, North: -50, Up: 20},
		"B": {East: 110, North: -48, Up: 20.5},
		"C": {East: 108, North: -38, Up: 21},
		"D": {East: 98, North: -40, Up: 19.5},
	}
}

func networkProblem(name string, truth map[string]model.ENU, fixedID string) model.Problem {
	order := []string{"A", "B", "C", "D", "E"}
	problem := model.Problem{Name: name}
	for _, id := range order {
		position, exists := truth[id]
		if !exists {
			continue
		}
		station := model.Station{ID: id}
		if id == fixedID {
			known := position
			station.Fixed = true
			station.KnownENU = &known
		}
		problem.Stations = append(problem.Stations, station)
	}
	return problem
}

func addBaselines(problem *model.Problem, truth map[string]model.ENU, fixtures []baselineFixture, covariance model.Matrix3) {
	for _, fixture := range fixtures {
		vector := addENU(subtractENU(truth[fixture.To], truth[fixture.From]), fixture.Noise)
		problem.Baselines = append(problem.Baselines, model.Baseline{
			ID: fixture.ID, From: fixture.From, To: fixture.To, Vector: vector,
			Covariance: covariance, Group: fixture.Group,
		})
	}
}

func patternedNoise(index int, scale float64) model.ENU {
	return model.ENU{
		East:  scale * float64(index%5-2),
		North: scale * float64((index*2)%7-3),
		Up:    scale * 0.5 * float64((index*3)%5-2),
	}
}

func addENU(left, right model.ENU) model.ENU {
	return model.ENU{East: left.East + right.East, North: left.North + right.North, Up: left.Up + right.Up}
}

func subtractENU(left, right model.ENU) model.ENU {
	return model.ENU{East: left.East - right.East, North: left.North - right.North, Up: left.Up - right.Up}
}
