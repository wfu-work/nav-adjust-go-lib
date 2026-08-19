package sparse

import (
	"fmt"
	"math"
)

// Preconditioner approximately solves M*x=rhs for one PCG iteration.
// Implementations must permit destination and rhs to refer to different
// slices of the matrix dimension.
type Preconditioner interface {
	Apply(destination, rhs []float64) error
	Name() string
}

type jacobiPreconditioner struct {
	inverseDiagonal []float64
}

// NewJacobiPreconditioner builds the scalar diagonal preconditioner used by
// the original sparse solver. Project may identify zero diagonal parameters
// that are completely removed by exact constraints.
func NewJacobiPreconditioner(matrix *Matrix, project ProjectFunc) (Preconditioner, error) {
	if matrix == nil {
		return nil, fmt.Errorf("sparse: Jacobi preconditioner requires a matrix")
	}
	inverse := make([]float64, matrix.Size())
	for i, value := range matrix.diagonal {
		if value < 0 || !finite(value) {
			return nil, fmt.Errorf("%w: invalid diagonal at %d", ErrBreakdown, i)
		}
		if value == 0 {
			if err := validateProjectedZero(matrix.Size(), i, project); err != nil {
				return nil, err
			}
			inverse[i] = 1
			continue
		}
		inverse[i] = 1 / value
	}
	return &jacobiPreconditioner{inverseDiagonal: inverse}, nil
}

func (preconditioner *jacobiPreconditioner) Apply(destination, rhs []float64) error {
	if preconditioner == nil || len(destination) != len(preconditioner.inverseDiagonal) || len(rhs) != len(destination) {
		return fmt.Errorf("sparse: Jacobi preconditioner dimension mismatch")
	}
	for i, value := range rhs {
		destination[i] = preconditioner.inverseDiagonal[i] * value
	}
	return nil
}

func (*jacobiPreconditioner) Name() string { return "jacobi" }

type diagonalBlock struct {
	start  int
	size   int
	factor []float64
	zero   bool
}

type blockJacobiPreconditioner struct {
	size   int
	blocks []diagonalBlock
}

// NewBlockJacobiPreconditioner factorizes consecutive principal diagonal
// blocks. For an ENU network blockSize=3 groups each station's East, North,
// and Up parameters while any final partial block remains valid.
func NewBlockJacobiPreconditioner(matrix *Matrix, blockSize int, project ProjectFunc) (Preconditioner, error) {
	if matrix == nil {
		return nil, fmt.Errorf("sparse: block-Jacobi preconditioner requires a matrix")
	}
	if blockSize <= 0 {
		return nil, fmt.Errorf("sparse: block-Jacobi block size must be positive")
	}
	result := &blockJacobiPreconditioner{size: matrix.Size()}
	for start := 0; start < matrix.Size(); start += blockSize {
		size := min(blockSize, matrix.Size()-start)
		factor := make([]float64, size*size)
		allZero := true
		for row := range size {
			for column := 0; column <= row; column++ {
				value := matrix.At(start+row, start+column)
				if !finite(value) {
					return nil, fmt.Errorf("%w: non-finite diagonal block at %d", ErrBreakdown, start)
				}
				factor[row*size+column] = value
				factor[column*size+row] = value
				allZero = allZero && value == 0
			}
		}
		block := diagonalBlock{start: start, size: size, factor: factor}
		if allZero {
			for index := range size {
				if err := validateProjectedZero(matrix.Size(), start+index, project); err != nil {
					return nil, err
				}
			}
			block.zero = true
			result.blocks = append(result.blocks, block)
			continue
		}
		if err := choleskyFactor(block.factor, size); err != nil {
			return nil, fmt.Errorf("%w: diagonal block starting at %d is not positive definite", ErrBreakdown, start)
		}
		result.blocks = append(result.blocks, block)
	}
	return result, nil
}

func (preconditioner *blockJacobiPreconditioner) Apply(destination, rhs []float64) error {
	if preconditioner == nil || len(destination) != preconditioner.size || len(rhs) != preconditioner.size {
		return fmt.Errorf("sparse: block-Jacobi preconditioner dimension mismatch")
	}
	for _, block := range preconditioner.blocks {
		if block.zero {
			copy(destination[block.start:block.start+block.size], rhs[block.start:block.start+block.size])
			continue
		}
		solveCholesky(destination[block.start:block.start+block.size], rhs[block.start:block.start+block.size], block.factor, block.size)
	}
	return nil
}

func (*blockJacobiPreconditioner) Name() string { return "block-jacobi" }

type incompleteCholeskyPreconditioner struct {
	size            int
	rowStart        []int
	columns         []int
	values          []float64
	diagonalOffsets []int
}

// NewIncompleteCholeskyPreconditioner builds an IC(0) factor using the lower
// triangle's existing non-zero pattern. relativeShift adds a small diagonal
// stabilization relative to each original diagonal entry; it changes only the
// preconditioner, not the system solved by PCG.
func NewIncompleteCholeskyPreconditioner(matrix *Matrix, relativeShift float64) (Preconditioner, error) {
	if matrix == nil {
		return nil, fmt.Errorf("sparse: IC(0) preconditioner requires a matrix")
	}
	if relativeShift < 0 || !finite(relativeShift) {
		return nil, fmt.Errorf("sparse: IC(0) relative shift must be finite and non-negative")
	}

	result := &incompleteCholeskyPreconditioner{
		size:            matrix.Size(),
		rowStart:        make([]int, matrix.Size()+1),
		columns:         make([]int, 0, len(matrix.columns)/2+matrix.Size()),
		values:          make([]float64, 0, len(matrix.values)/2+matrix.Size()),
		diagonalOffsets: make([]int, matrix.Size()),
	}
	for row := 0; row < matrix.Size(); row++ {
		foundDiagonal := false
		for offset := matrix.rowStart[row]; offset < matrix.rowStart[row+1]; offset++ {
			column := matrix.columns[offset]
			if column > row {
				break
			}
			result.columns = append(result.columns, column)
			result.values = append(result.values, matrix.values[offset])
			if column == row {
				foundDiagonal = true
			}
		}
		if !foundDiagonal {
			result.columns = append(result.columns, row)
			result.values = append(result.values, 0)
		}
		result.diagonalOffsets[row] = len(result.values) - 1
		result.rowStart[row+1] = len(result.values)
	}

	for row := 0; row < result.size; row++ {
		rowStart := result.rowStart[row]
		diagonalOffset := result.diagonalOffsets[row]
		for offset := rowStart; offset < diagonalOffset; offset++ {
			column := result.columns[offset]
			value := result.values[offset]
			left := rowStart
			right := result.rowStart[column]
			rightEnd := result.diagonalOffsets[column]
			for left < offset && right < rightEnd {
				leftColumn := result.columns[left]
				rightColumn := result.columns[right]
				switch {
				case leftColumn < rightColumn:
					left++
				case leftColumn > rightColumn:
					right++
				default:
					value -= result.values[left] * result.values[right]
					left++
					right++
				}
			}
			pivot := result.values[result.diagonalOffsets[column]]
			if pivot <= 0 || !finite(pivot) {
				return nil, fmt.Errorf("%w: IC(0) invalid pivot at %d", ErrBreakdown, column)
			}
			result.values[offset] = value / pivot
			if !finite(result.values[offset]) {
				return nil, fmt.Errorf("%w: IC(0) non-finite factor at [%d,%d]", ErrBreakdown, row, column)
			}
		}

		originalDiagonal := result.values[diagonalOffset]
		pivot := originalDiagonal + relativeShift*math.Max(math.Abs(originalDiagonal), 1)
		for offset := rowStart; offset < diagonalOffset; offset++ {
			pivot -= result.values[offset] * result.values[offset]
		}
		if pivot <= 0 || !finite(pivot) {
			return nil, fmt.Errorf("%w: IC(0) non-positive pivot at %d", ErrBreakdown, row)
		}
		result.values[diagonalOffset] = math.Sqrt(pivot)
	}
	return result, nil
}

func (preconditioner *incompleteCholeskyPreconditioner) Apply(destination, rhs []float64) error {
	if preconditioner == nil || len(destination) != preconditioner.size || len(rhs) != preconditioner.size {
		return fmt.Errorf("sparse: IC(0) preconditioner dimension mismatch")
	}
	for row := 0; row < preconditioner.size; row++ {
		value := rhs[row]
		for offset := preconditioner.rowStart[row]; offset < preconditioner.diagonalOffsets[row]; offset++ {
			value -= preconditioner.values[offset] * destination[preconditioner.columns[offset]]
		}
		destination[row] = value / preconditioner.values[preconditioner.diagonalOffsets[row]]
	}
	for row := preconditioner.size - 1; row >= 0; row-- {
		diagonal := preconditioner.values[preconditioner.diagonalOffsets[row]]
		destination[row] /= diagonal
		for offset := preconditioner.rowStart[row]; offset < preconditioner.diagonalOffsets[row]; offset++ {
			column := preconditioner.columns[offset]
			destination[column] -= preconditioner.values[offset] * destination[row]
		}
	}
	return nil
}

func (*incompleteCholeskyPreconditioner) Name() string { return "ic0" }

func validateProjectedZero(size, index int, project ProjectFunc) error {
	if project == nil {
		return fmt.Errorf("%w: invalid diagonal at %d", ErrBreakdown, index)
	}
	basis := make([]float64, size)
	basis[index] = 1
	if err := project(basis); err != nil {
		return fmt.Errorf("sparse: test projected diagonal %d: %w", index, err)
	}
	if norm(basis) > 1e-12 {
		return fmt.Errorf("%w: unconstrained zero diagonal at %d", ErrBreakdown, index)
	}
	return nil
}

func choleskyFactor(matrix []float64, size int) error {
	for row := range size {
		for column := 0; column <= row; column++ {
			value := matrix[row*size+column]
			for k := 0; k < column; k++ {
				value -= matrix[row*size+k] * matrix[column*size+k]
			}
			if row == column {
				if value <= 0 || !finite(value) {
					return ErrBreakdown
				}
				matrix[row*size+column] = math.Sqrt(value)
			} else {
				matrix[row*size+column] = value / matrix[column*size+column]
			}
		}
		for column := row + 1; column < size; column++ {
			matrix[row*size+column] = 0
		}
	}
	return nil
}

func solveCholesky(destination, rhs, factor []float64, size int) {
	for row := range size {
		value := rhs[row]
		for column := 0; column < row; column++ {
			value -= factor[row*size+column] * destination[column]
		}
		destination[row] = value / factor[row*size+row]
	}
	for row := size - 1; row >= 0; row-- {
		value := destination[row]
		for column := row + 1; column < size; column++ {
			value -= factor[column*size+row] * destination[column]
		}
		destination[row] = value / factor[row*size+row]
	}
}
