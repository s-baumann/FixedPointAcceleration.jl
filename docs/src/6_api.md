```@meta
CurrentModule = FixedPointAcceleration
```

# Function API

```@index
Pages = ["api.md"]
```

### Main fixed point seeking function

```@docs
    fixed_point
    fixed_point_new_input
```

### Structs

```@docs
    FunctionEvaluationResult
    FixedPointResults
```

### Internal Functions

```@docs
    put_together_without_jumps
    execute_function_safely
```

### Algorithm-Specific Internal Functions

The extrapolation methods are now implemented directly within each algorithm:
- **MPE**: `_mpe_extrapolation()` in `src/algorithms/mpe.jl`
- **RRE**: `_rre_extrapolation()` in `src/algorithms/rre.jl` 
- **VEA**: `_vea_epsilon_extrapolation()` and `_vea_vector_of_inverses()` in `src/algorithms/vea.jl`
- **SEA**: `_sea_epsilon_extrapolation()` in `src/algorithms/sea.jl`
