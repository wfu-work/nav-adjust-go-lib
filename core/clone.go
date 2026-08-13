package core

// CloneProblem returns a deep copy suitable for iterative reweighting or model
// adaptation without mutating the caller's input.
func CloneProblem(src Problem) Problem {
	dst := src
	dst.Parameters = append([]Parameter(nil), src.Parameters...)
	dst.Equations = append([]Equation(nil), src.Equations...)
	for i := range dst.Equations {
		dst.Equations[i].Terms = append([]Term(nil), src.Equations[i].Terms...)
	}
	dst.CovarianceBlocks = append([]CovarianceBlock(nil), src.CovarianceBlocks...)
	for i := range dst.CovarianceBlocks {
		dst.CovarianceBlocks[i].RowIndexes = append([]int(nil), src.CovarianceBlocks[i].RowIndexes...)
		dst.CovarianceBlocks[i].Covariance = append([]float64(nil), src.CovarianceBlocks[i].Covariance...)
	}
	dst.Constraints = append([]Constraint(nil), src.Constraints...)
	for i := range dst.Constraints {
		dst.Constraints[i].Terms = append([]Term(nil), src.Constraints[i].Terms...)
	}
	return dst
}
