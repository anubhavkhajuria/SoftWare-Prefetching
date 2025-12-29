# Software Prefetching for Irregular Memory Access

An LLVM pass that automatically inserts prefetch instructions for irregular memory access patterns like `A[B[i]]`, commonly found in graph algorithms.

## Features

- **Automatic IMA Detection**: Identifies indirect array accesses through static slicing analysis
- **Cost Model**: Scores loops to skip unprofitable prefetches
- **Graph Reordering**: Includes Gorder, Hub Sorting, and Degree-Based Grouping algorithms
- **Configurable**: Tunable prefetch distance, locality hints, and batch size

## Quick Start

```bash
# Build the LLVM pass
./build.sh

# Run PageRank benchmarks with graph reordering
./run_tests.sh rcm      # Reverse Cuthill-McKee
./run_tests.sh gorder   # Gorder algorithm
./run_tests.sh hub      # Hub Sorting
./run_tests.sh none     # No reordering
```


## Pass Options

| Option | Default | Description |
|--------|---------|-------------|
| `--ima-prefetch-distance` | 8 | Iterations to prefetch ahead |
| `--ima-degree-threshold` | 4 | Minimum loop iterations |
| `--ima-locality` | 3 | Cache hint (3=L1, 2=L2, 1=L3) |
| `--ima-batch-size` | 1 | Addresses per prefetch (1-8) |
| `--ima-cost-threshold` | 20 | Minimum benefit score (0-100) |
| `--ima-enable-cost-model` | true | Enable cost-based filtering |

## Usage

```bash
# Compile to LLVM IR
clang++ -emit-llvm -S input.cpp -o input.ll -std=c++17

# Preprocess
opt -passes='mem2reg,loop-simplify,lcssa,indvars' -S input.ll -o input_prep.ll

# Apply prefetch pass
opt -load-pass-plugin ./build/pass/prefetchIMAPass.dylib \
    -passes='prefetch-ima' \
    -ima-prefetch-distance=8 \
    input_prep.ll -o output.ll

# Compile final executable
clang++ -O3 output.ll -o output
```

## Graph Reordering Algorithms

The testcases support multiple reordering strategies in `gorder.h`:

- **RCM**: Reverse Cuthill-McKee (bandwidth reduction)
- **Gorder**: Greedy window-based neighborhood overlap
- **Hub Sorting**: Separates high-degree hubs from cold vertices
- **DBG**: Degree-Based Grouping into buckets
- **Hub+Gorder**: Combined approach

## Requirements

- LLVM 14+ with development headers
- CMake 3.12+
- C++17 compiler

## Platform Notes

- **macOS**: Uses `prefetchIMAPass.dylib`
- **Linux**: Uses `prefetchIMAPass.so`

Update `run_tests.sh` accordingly for your platform.
