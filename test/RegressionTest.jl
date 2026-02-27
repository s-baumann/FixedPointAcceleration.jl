# Regression tests for FixedPointAcceleration.jl
#
# Reference values live in regression_refs.yaml, which is also consumed by the
# Rust (and any future language) implementation for cross-language comparison.
#
# ── Normal test run ───────────────────────────────────────────────────────────
#   julia --project -e 'using Pkg; Pkg.test()'
#
# ── Rebase (update expected values after an intentional algorithm change) ─────
#   julia -e 'using Pkg; withenv("REBASE_REFS" => "true") do; Pkg.test(); end'
#
# On rebase the following fields are overwritten per algorithm:
#   expected_termination, expected_fixed_point, iterates (first 4 input/output pairs)
# expected_iterations is deleted if still present from an older schema.
# iterate_tol is added at test-case level if missing (default 1e-12).
#
# All other fields (function, x0, convergence_threshold, fixed_point_tol, …)
# are human-authored and are never touched by rebase.

using YAML, Test, FixedPointAcceleration

const REFS_FILE = joinpath(@__DIR__, "regression_refs.yaml")
const REBASE    = get(ENV, "REBASE_REFS", "false") == "true"

# ── Function registry ─────────────────────────────────────────────────────────
# Functions are identified by the string key in the YAML `function` field.
# Every language that consumes regression_refs.yaml must maintain a matching
# registry with the same keys and equivalent implementations.
const FUNC_REGISTRY = Dict{String, Function}(
    "cos"               => x -> cos.(x),
    "vector_2d"         => x -> [0.5 * sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]],
    "complex_linear_1d" => z -> (z .+ ComplexF64(1, 2)) ./ 2,
    "complex_linear_2d" => v -> (v .+ ComplexF64[1+2im, -3+1im]) ./ 2,
    "sqrt"              => x -> sqrt.(x),
)

# ── YAML ↔ Julia helpers ──────────────────────────────────────────────────────

# Decode a raw YAML value (list or {real,imag} dict) into a Julia vector.
function decode_fp(val)
    val isa Dict ? ComplexF64.(val["real"], val["imag"]) : Float64.(val)
end

# Build x0 from the YAML x0 spec.
build_x0(spec) = decode_fp(spec["x0"])

# Encode a fixed-point or iterate vector for writing back to YAML.
function encode_fp(fp::AbstractVector)
    eltype(fp) <: Complex ?
        Dict("real" => collect(Float64, real.(fp)), "imag" => collect(Float64, imag.(fp))) :
        collect(Float64, fp)
end

# Build keyword arguments from the test-case spec.
function build_kwargs(spec)
    opts = get(spec, "options", nothing)
    kw = (
        MaxIter                    = Int(spec["max_iter"]),
        ConvergenceMetricThreshold = Float64(spec["convergence_threshold"]),
    )
    if !isnothing(opts)
        kw = merge(kw, (
            Dampening            = Float64(get(opts, "dampening", 1.0)),
            Dampening_With_Input = Bool(get(opts, "dampening_with_input", false)),
        ))
    end
    kw
end

# ── Core: run (or rebase) one test case ──────────────────────────────────────

function run_case!(spec)
    func     = FUNC_REGISTRY[spec["function"]]
    x0       = build_x0(spec)
    fp_tol   = Float64(spec["fixed_point_tol"])
    iter_tol = Float64(get(spec, "iterate_tol", 1e-12))
    conv_thr = Float64(spec["convergence_threshold"])
    kw       = build_kwargs(spec)
    algs     = spec["algorithms"]

    for (alg_name, alg_spec) in algs
        alg = Symbol(alg_name)
        r   = fixed_point(func, x0; Algorithm = alg, kw...)

        if REBASE
            # Ensure iterate_tol is present at the test-case level.
            if !haskey(spec, "iterate_tol")
                spec["iterate_tol"] = 1e-12
            end
            # Overwrite the three machine-generated fields.
            alg_spec["expected_termination"] = string(r.TerminationCondition_)
            alg_spec["expected_fixed_point"] = encode_fp(r.FixedPoint_)
            n_store = min(4, size(r.Inputs_, 2))
            alg_spec["iterates"] = [
                Dict("input"  => encode_fp(r.Inputs_[:,i]),
                     "output" => encode_fp(r.Outputs_[:,i]))
                for i in 1:n_store
            ]
            # Remove legacy field if present from an older schema version.
            delete!(alg_spec, "expected_iterations")
        else
            @testset "$alg_name" begin
                @test string(r.TerminationCondition_) == alg_spec["expected_termination"]
                @test r.FixedPoint_ ≈ decode_fp(alg_spec["expected_fixed_point"]) atol=fp_tol
                @test r.Convergence_ < conv_thr

                for (i, iter_ref) in enumerate(get(alg_spec, "iterates", []))
                    @testset "iterate $i" begin
                        @test r.Inputs_[:,i]  ≈ decode_fp(iter_ref["input"])  atol=iter_tol
                        @test r.Outputs_[:,i] ≈ decode_fp(iter_ref["output"]) atol=iter_tol
                    end
                end
            end
        end
    end
end

# ── Entry point ───────────────────────────────────────────────────────────────

refs = YAML.load_file(REFS_FILE)

if REBASE
    for spec in refs["test_cases"]
        run_case!(spec)
    end
    open(REFS_FILE, "w") do io
        YAML.write(io, refs)
    end
    @info "Rebased → $REFS_FILE  (commit the updated file)"
else
    @testset "Regression Tests" begin
        for spec in refs["test_cases"]
            @testset "$(spec["name"])" begin
                run_case!(spec)
            end
        end
    end
end
