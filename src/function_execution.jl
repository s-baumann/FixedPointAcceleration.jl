"""
    execute_function_safely(Func::Function, x::Array{T,1}; type_check::Bool = false, quiet_errors::Bool = true) where T<:Real
This function creates a function that executes the function for which a fixed point is sought. It is a helper function that is not exported.
### Inputs
 * `Func` - The function input to fixed_point
 * `x` - The point at which to evaluate the function.
 * `type_check` - Should the type stability of the function be checked and a error reported in an input vector of Ints turns into a vector of floats (for instance)
### Returns
 * A `FunctionEvaluationResult`
### Examples
    Func(x) = sqrt(x)
    execute_function_safely(Func, [-1.0,0.0,1.0])
    execute_function_safely(Func,[Missing(),0.0,1.0])
    execute_function_safely(Func,[7.0,0.0,1.0])
    execute_function_safely(Func,[NaN,0.0,1.0])
    execute_function_safely(Func,[Inf,0.0,1.0])
    execute_function_safely(Func,-1.0)
    execute_function_safely(Func,Missing())
    execute_function_safely(Func,1.0)
    execute_function_safely(Func,NaN)
    execute_function_safely(Func,Inf)
"""
function execute_function_safely(Func::Function, x::Array{T,1}; type_check::Bool = false, quiet_errors::Bool = true) where T<:Real
    # Check input
    if sum(isnan.(x)) > 0
        return FunctionEvaluationResult(x, missing, :InputNAsDetected)
    elseif sum(isinf.(x)) > 0
        return FunctionEvaluationResult(x, missing, :InputInfsDetected)
    end
    # Run function.
    lenx = length(x)
    tf_result = Array{T,1}(undef,lenx)
    tf_full_result = missing
    if quiet_errors
        try
            tf_full_result =  Func(deepcopy(x))
        catch
            return FunctionEvaluationResult(x, missing, :ErrorExecutingFunction)
        end
    else
        tf_full_result =  Func(deepcopy(x))
    end

    side_effect_to_report = missing
    # Now we if tf_full_result is an array then it is the normal case.
    if isa(tf_full_result,Vector)
        tf_result = tf_full_result
    elseif isa(tf_full_result,Tuple) && (length(tf_full_result) == 2) && ((isa(tf_full_result[1],Vector) & isa(tf_full_result[2],NamedTuple)))
        tf_result = tf_full_result[1]
        side_effect_to_report = tf_full_result[2]
    else
        error("This function returned a $(typeof(tf_full_result)). The Fixedpoint function can only return a vector or a tuple of which the first entry is the vector for which a fixedpoint is sought and the second is a namedtuple (the contents of which are output for the user but are not used in fixed point acceleration).")
    end
    # Check Output and return.
    if sum(ismissing.(tf_result)) > 0
        return FunctionEvaluationResult(x, tf_result, :OutputMissingsDetected      , side_effect_to_report)
    elseif sum(isnan.(tf_result)) > 0
        return FunctionEvaluationResult(x, tf_result, :OutputNAsDetected           , side_effect_to_report)
    elseif sum(isinf.(tf_result)) > 0
        return FunctionEvaluationResult(x, tf_result, :OutputInfsDetected          , side_effect_to_report)
    elseif (length(tf_result) != length(x))
        return FunctionEvaluationResult(x, tf_result, :LengthOfOutputNotSameAsInput, side_effect_to_report)
    elseif type_check & (!isa(tf_result[1], T))
        return FunctionEvaluationResult(x, tf_result, :FunctionIsNotTypeStable     , side_effect_to_report)
    else
        return FunctionEvaluationResult(x, tf_result, :NoError                     , side_effect_to_report)
    end
end
