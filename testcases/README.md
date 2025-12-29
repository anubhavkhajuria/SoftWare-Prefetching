# Testcases and Preprocessing

This directory contains PageRank implementations and preprocessing utilities for the PrefetchIMA pass.

## Files

| File | Description |
|------|-------------|
| `pagerank_array.cpp` | PageRank using C-style arrays (faster) |
| `pagerank_vector.cpp` | PageRank using std::vector |
| `gorder.h` | Graph reordering algorithms |
| `prefetch_stats.h` | Runtime profiling utilities |

## Graph Reordering

Both PageRank implementations support 6 reordering modes:

```bash
./pagerank_array_exe <graph_file> <start_node> <max_nodes> [mode]
```

| Mode | Description |
|------|-------------|
| `none` | No reordering (baseline) |
| `rcm` | Reverse Cuthill-McKee |
| `gorder` | Greedy window overlap |
| `hub` | Hub Sorting (default) |
| `dbg` | Degree-Based Grouping |
| `hubgorder` | Hub Sorting + Gorder |

## Preprocessing Pipeline

The compilation pipeline applies LLVM passes before the prefetch pass:

```bash
# 1. Emit LLVM IR
clang++ -emit-llvm -S input.cpp -o input.ll -std=c++17 -O0 -fno-inline

# 2. Apply canonicalization passes
opt -passes='mem2reg,loop-simplify,lcssa,indvars,simplifycfg' \
    -S input.ll -o input_prep.ll

# 3. Apply prefetch pass
opt -load-pass-plugin ./build/pass/prefetchIMAPass.dylib \
    -passes='prefetch-ima' \
    -S input_prep.ll -o input_prefetch.ll

# 4. Compile to optimized binary
clang++ -O3 input_prefetch.ll -o output
```

### Required LLVM Passes

| Pass | Purpose |
|------|---------|
| `mem2reg` | Promote memory to registers (SSA form) |
| `loop-simplify` | Canonicalize loop structure |
| `lcssa` | Loop-Closed SSA form |
| `indvars` | Canonicalize induction variables |
| `simplifycfg` | Simplify control flow |

## Runtime Profiling

Enable profiling by defining `PROFILE_PREFETCH`:

```cpp
#define PROFILE_PREFETCH
#include "prefetch_stats.h"

int main() {
    // ... your code ...
    PREFETCH_PRINT_STATS();
    return 0;
}
```

Output includes prefetch count, overhead timing, and memory statistics.

## Input Format

Graph files use edge list format (SNAP):

```
# Optional comments starting with #
<from_node> <to_node>
<from_node> <to_node>
...
```

Example datasets available from [SNAP](https://snap.stanford.edu/data/).
