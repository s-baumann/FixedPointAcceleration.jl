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
	EpsilonExtrapolationVectorOfInverses
	PolynomialExtrapolation
	EpsilonExtrapolation
```
