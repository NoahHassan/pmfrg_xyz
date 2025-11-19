# Performance Engineering for pmfrg_xyz

This Julia project contains tools for profiling and benchmarking the code in pmfrg_xyz.

## Installed Packages

- **PProf.jl**: Statistical profiling with pprof integration
- **BenchmarkTools.jl**: Accurate benchmarking framework
- **ProfileSVG.jl**: SVG-based profile visualization
- **FlameGraphs.jl**: Flame graph visualization for profiling

## Usage
Activate this project environment:
```bash
julia --project=pmfrg_xyz/performance_engineering
```

### Profiling with PProf

```julia
using PProf

# Profile your code
@pprof myfunction(args...)

# Or profile a code block
@profile begin
    # your code here
end
pprof()
```

### Benchmarking with BenchmarkTools

```julia
using BenchmarkTools

# Quick benchmark
@btime myfunction($args)

# Detailed benchmark
@benchmark myfunction($args)

# Save benchmark results
b = @benchmark myfunction($args)
```

### Profile Visualization

```julia
using ProfileSVG, Profile

# Profile and save as SVG
@profile myfunction(args...)
ProfileSVG.save("profile.svg")

# Or use FlameGraphs
using FlameGraphs
@profile myfunction(args...)
flamegraph()
```

## Ready-to-Use Scripts

### Profile Dimer Anisotropy Example

Profile the dimer anisotropy calculation with PProf:
```bash
julia pmfrg_xyz/performance_engineering/profile_dimer_anisotropy.jl
```

This will:
- Run the calculation once for compilation
- Profile the second run with detailed call stack information
- Open an interactive pprof viewer in your browser

### Benchmark Dimer Anisotropy Example

Benchmark the dimer anisotropy calculation:
```bash
julia pmfrg_xyz/performance_engineering/benchmark_dimer_anisotropy.jl
```

This will:
- Run the calculation once for compilation
- Perform multiple runs to get accurate timing statistics
- Display min/median/mean/max execution times

## Example Workflow

1. Use the ready-to-use scripts above for common profiling tasks
2. Or create custom profiling scripts:
   - Activate the environment: `julia --project=pmfrg_xyz/performance_engineering`
   - Load your code: `include("../Lpmfrg_xyz.jl")` or `include("../Tpmfrg_xyz.jl")`
   - Profile or benchmark the functions you want to optimize
   - Visualize results and identify bottlenecks
