// Package statistics provides internal post-adjustment statistical helpers.
package statistics

import (
	"fmt"
	"math"
	"sort"

	adjust "github.com/wfu-work/nav-adjust-go-lib/core"
	"gonum.org/v1/gonum/stat/distuv"
)

// GlobalTest describes a two-sided chi-square test of the a-priori variance
// factor (expected value 1).
type GlobalTest struct {
	Statistic  float64
	DOF        int
	Confidence float64
	Lower      float64
	Upper      float64
	Passed     bool
}

// ChiSquare tests Result.Objective against a chi-square distribution.
func ChiSquare(result *adjust.Result, confidence float64) (GlobalTest, error) {
	if result == nil || result.DOF <= 0 {
		return GlobalTest{}, fmt.Errorf("statistics: chi-square test requires positive degrees of freedom")
	}
	if result.Objective < 0 || math.IsNaN(result.Objective) || math.IsInf(result.Objective, 0) {
		return GlobalTest{}, fmt.Errorf("statistics: chi-square test requires a finite non-negative objective")
	}
	if confidence <= 0 || confidence >= 1 || math.IsNaN(confidence) || math.IsInf(confidence, 0) {
		return GlobalTest{}, fmt.Errorf("statistics: confidence must be between zero and one")
	}
	distribution := distuv.ChiSquared{K: float64(result.DOF)}
	alpha := 1 - confidence
	lower := distribution.Quantile(alpha / 2)
	upper := distribution.Quantile(1 - alpha/2)
	return GlobalTest{
		Statistic: result.Objective, DOF: result.DOF, Confidence: confidence,
		Lower: lower, Upper: upper,
		Passed: result.Objective >= lower && result.Objective <= upper,
	}, nil
}

// ResidualCandidate is an observation sorted by absolute standardized residual.
type ResidualCandidate struct {
	ID           string
	Group        string
	Standardized float64
}

// LargestResiduals returns at most limit residuals in descending absolute
// standardized-residual order. A non-positive limit returns all residuals.
func LargestResiduals(result *adjust.Result, limit int) []ResidualCandidate {
	if result == nil {
		return nil
	}
	candidates := make([]ResidualCandidate, len(result.Residuals))
	for i, residual := range result.Residuals {
		candidates[i] = ResidualCandidate{residual.ID, residual.Group, residual.Standardized}
	}
	sort.SliceStable(candidates, func(i, j int) bool {
		return math.Abs(candidates[i].Standardized) > math.Abs(candidates[j].Standardized)
	})
	if limit > 0 && limit < len(candidates) {
		candidates = candidates[:limit]
	}
	return candidates
}
