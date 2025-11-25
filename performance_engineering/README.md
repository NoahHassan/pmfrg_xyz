# Performance Engineering for pmfrg_xyz

This Julia project contains tools for profiling and benchmarking the code in pmfrg_xyz.

## Installed Packages

- **PProf.jl**: Statistical profiling with pprof integration
- **BenchmarkTools.jl**: Accurate benchmarking framework
- **ProfileSVG.jl**: SVG-based profile visualization
- **FlameGraphs.jl**: Flame graph visualization for profiling

## Ready-to-Use Scripts

### Profiling Examples

Profile any example with PProf (`dimer` or `square_lattice`):
```bash
julia pmfrg_xyz/performance_engineering/profile.jl [example_name]
```

Examples:
```bash
# Profile dimer example (default)
julia pmfrg_xyz/performance_engineering/profile.jl
julia pmfrg_xyz/performance_engineering/profile.jl dimer

# Profile square lattice example
julia pmfrg_xyz/performance_engineering/profile.jl square_lattice
```

This will:
- Run the calculation once for compilation
- Profile the second run with detailed call stack information
- Save profile data to `profile_data/profile_<example>_<commit>.pb.gz`
- Open an interactive pprof viewer in your browser

### Benchmarking Examples

#### Full Example Benchmarks (WIP)

The script `pmfrg_xyz/performance_engineering/benchmark.jl`
works in similar way as the profiling script and:
- Perform multiple runs to get accurate timing statistics
- Display min/median/mean/max execution times
- Save results to `benchmark_data/benchmark_<example>_<commit>.txt`

Note: this is now impractical because these examples
take still too long
for their benchmark to be useful.

#### getXBubble! Function Benchmarks

The script `benchmark_getXBubble.jl` provides focused benchmarks for the `getXBubble!` function:

```bash
julia --project=performance_engineering performance_engineering/benchmark_getXBubble.jl
```

This script runs three types of benchmarks:

1. **Regression test data**: Uses real data from the test suite
2. **Synthetic dimer data**: Parametrized by `N` (frequency discretization)
   - Default: `N=8`
   - Modify in script to test different values
3. **Synthetic square lattice data**: Parametrized by `N` and `lattice_size`
   - Default: `N=10`, `lattice_size=8`
   - Modify in script to test different system sizes

Each benchmark provides:
- Timing distribution (min/median/mean/max)
- Memory usage and allocation counts
- Visual histogram of execution times

**Creating synthetic benchmarks**: The script includes helper functions to generate test data:
- `create_synthetic_workspace_dimer(N)`: Creates dimer system workspace
- `create_synthetic_workspace_square(N, lattice_size)`: Creates square lattice workspace 

