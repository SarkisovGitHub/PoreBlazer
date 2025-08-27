# OpenMP Parallelization in PoreBlazer 4.0

## Overview

This document describes the OpenMP parallelization implemented in PoreBlazer 4.0 to improve performance on multi-core systems.

## Parallelized Components

### 1. Lattice Calculations (`lattice_calculations` subroutine)

**Location**: `src/poreblazer.f90`, lines ~480-580

**Description**: The most computationally intensive part of the code, involving a triple nested loop over lattice cubes and an inner loop over all atoms.

**Parallelization Strategy**:
- Uses `!$OMP PARALLEL DO` with `COLLAPSE(3)` to parallelize the three nested loops over lattice dimensions
- Each thread processes different lattice cubes independently
- Dynamic scheduling for load balancing
- Thread-safe collection of results using sequential post-processing

**Performance Impact**: This is where the largest speedup is expected since it's the computational bottleneck.

### 2. Surface Area Calculations (`surface_area` subroutine)

**Location**: `src/poreblazer.f90`, lines ~815-910

**Description**: Monte Carlo sampling for accessible surface area calculation, involving loops over atoms and samples.

**Parallelization Strategy**:
- Uses `!$OMP PARALLEL DO` to parallelize the outer loop over atoms
- Each thread processes different atoms independently
- Critical sections protect random number generation to ensure thread safety
- Reduction clause for accumulating total surface area (`stotal`)

**Performance Impact**: Good parallelization potential with minimal synchronization overhead.

### 3. Data Collection Loops

**Locations**: 
- Limiting diameter analysis: lines ~1210-1230
- Nitrogen accessible cubes: lines ~1015-1035

**Description**: Loops that collect and organize data for further analysis.

**Parallelization Strategy**:
- Uses `!$OMP PARALLEL DO` with `COLLAPSE(3)` for nested loops
- Critical sections protect shared counter variables and array assignments
- Dynamic scheduling for load balancing

## Compilation Requirements

The code requires OpenMP support in the compiler. The Makefile has been updated to include:

```makefile
F90FLAGS = -fopenmp
OFLAGS = -O2 -unshared -fopenmp  
LINKERFLAGS = -fopenmp
```

## Usage

### Setting Number of Threads

Control the number of OpenMP threads using the `OMP_NUM_THREADS` environment variable:

```bash
export OMP_NUM_THREADS=4    # Use 4 threads
./poreblazer.exe input.dat
```

### Recommended Thread Count

- Start with the number of physical CPU cores
- For systems with hyperthreading, using all logical cores may or may not improve performance
- Test different thread counts to find optimal performance for your system

## Performance Considerations

### Thread Safety

- Random number generation is protected with critical sections
- Shared variables are properly handled with reduction clauses or critical sections
- Each thread works on independent data where possible

### Expected Speedup

- **Lattice calculations**: Near-linear speedup with thread count (limited by memory bandwidth)
- **Surface area calculations**: Good speedup, limited by random number generation synchronization
- **Data collection loops**: Moderate speedup due to critical section overhead

### Memory Requirements

- Parallel execution may slightly increase memory usage due to OpenMP overhead
- No significant additional memory allocation for parallelization

## Validation

The parallelized version has been tested to ensure:

1. **Correctness**: Core geometric properties (pore limiting diameter, maximum pore diameter, percolation dimensions) remain identical
2. **Consistency**: Small variations in Monte Carlo calculations (surface area) are expected and acceptable
3. **Stability**: No race conditions or memory corruption

## Notes

- The parallelization is conservative and safe, prioritizing correctness over maximum performance
- Random number generation could be further optimized with thread-local generators in future versions
- The current implementation maintains backward compatibility and identical algorithmic behavior

## Example Performance

Testing on a sample HKUST1 structure shows the parallelized version works correctly with both single and multiple threads, producing consistent results for deterministic calculations and acceptable variation for Monte Carlo components.