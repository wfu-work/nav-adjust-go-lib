package model

// PositionPrior is a stochastic control coordinate for one free station.
// Position and Covariance use the same common ENU frame as all baselines.
// Unlike Station.Fixed, a prior participates in the adjustment and may have a
// non-zero residual.
type PositionPrior struct {
	ID         string            `json:"id"`
	StationID  string            `json:"station_id"`
	Position   ENU               `json:"position"`
	Covariance Matrix3           `json:"covariance"`
	Metadata   map[string]string `json:"metadata,omitempty"`
}

// Prior is the short name used by applications that import model directly.
type Prior = PositionPrior
