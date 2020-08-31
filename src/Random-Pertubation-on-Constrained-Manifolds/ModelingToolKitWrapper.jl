using ModelingToolkit

"""
ModelingToolkit needs a special format to substitute values into operations.
This wrapper function inserts a new value into the Jacobi Matrix Jg.
Inputs:
        Jg: Symbolic Jacobi Matrix
        vars: variables of the Jacobi Matrix Jg
        eval_values: Values at which Jg should be evaluated
Outputs:
        Evaluated Jacobi Matrix
"""
function insert_value(Jg, vars, eval_values)
    num_vars = length(vars)
    # Inizilizes a Vector in the correct format for ModelingToolkit
    subsitute_vec = Vector{Pair{Operation,Float64}}(undef, num_vars)

    # Inserts eval_values as operations for vars
    for i in 1:num_vars
        subsitute_vec[i] = vars[i] => eval_values[i]
    end

    Jg_eval = substitute.(Jg,(subsitute_vec,)) # Evaluates Jg at eval_values
    Jg_eval = (p->p.value).(Jg_eval)  # Removes the ModelingToolkit.Constant

    Jg_eval = Float64.(Jg_eval) # Prevents that the Arry is of type Number
    return Jg_eval
end


"""
Calculates the symbolic Jacobian Jg of a function g using ModelingToolkit.
Takes a list of the arguments of g as symbols and turns them into Variables.
The variables and g are then  used to generate a experssion. The expression is
used to calculate the Jacobi Matrix Jg.
Inputs:
        g: Function we want the jacobi matrix from
        symbol_list: List containg the args of g as symbols
Outputs:
        Jg: Symbolic Jacobi Matrix
        vars: variables of the Jacobi Matrix Jg
"""
function jacobian_wrapper(g, symbol_list)
    num_vars = length(symbol_list)
    vars = Vector{Operation}(undef, num_vars)

    for k in 1:num_vars
        vars[k] = Variable(symbol_list[k])()
    end

    g_expr = g(vars...) # Using function g as a generator for an expression

    if (isa(g_expr ,Array{Operation})) == false
        g_expr = [g_expr] # ModelingToolkit needs an Array of Operations
    end

    Jg = ModelingToolkit.jacobian(g_expr, vars)
    Jg = simplify.(Jg)
    return Jg, vars
end
