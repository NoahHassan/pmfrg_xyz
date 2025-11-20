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

### Benchmarking Examples (STILL WIP)

The script `pmfrg_xyz/performance_engineering/benchmark.jl`
works in similar way as the profiling script and:
- Perform multiple runs to get accurate timing statistics
- Display min/median/mean/max execution times
- Save results to `benchmark_data/benchmark_<example>_<commit>.txt`

Note: this is now impractical because these examples 
take still too long 
for their benchmark to be useful.

TODO: Benchmarks would be useful for smaller functions. 

