package model

// Station describes a node in the ENU network. Fixed stations require
// KnownENU. Free stations are solved directly and do not require an initial
// position.
type Station struct {
	ID       string            `json:"id"`
	Name     string            `json:"name,omitempty"`
	Fixed    bool              `json:"fixed"`
	KnownENU *ENU              `json:"known_enu,omitempty"`
	Metadata map[string]string `json:"metadata,omitempty"`
}

// ENUBaseline is an observed relative vector From -> To. All baselines in a
// problem must use the same ENU axes and unit (metres).
type ENUBaseline struct {
	ID         string            `json:"id"`
	From       string            `json:"from"`
	To         string            `json:"to"`
	Vector     ENU               `json:"vector"`
	Covariance Matrix3           `json:"covariance"`
	Group      string            `json:"group,omitempty"`
	Metadata   map[string]string `json:"metadata,omitempty"`
}

// ENUNetworkProblem is the complete public ENU network input.
type ENUNetworkProblem struct {
	Name      string          `json:"name,omitempty"`
	Stations  []Station       `json:"stations"`
	Baselines []ENUBaseline   `json:"baselines"`
	Priors    []PositionPrior `json:"priors,omitempty"`
}

// Short input names are provided for applications that import model directly.
// The ENU-prefixed names remain the stable root-facade vocabulary.
type (
	Baseline = ENUBaseline
	Problem  = ENUNetworkProblem
)
