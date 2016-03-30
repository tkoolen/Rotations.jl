# defines common methods for the various rotation types
RotTypeList    = [RotMatrix, Quaternion, SpQuat, EulerAngles, ProperEulerAngles, AngleAxis]
RotationTypes  = Union{[RotMatrix, Quaternion, SpQuat, EulerAngles, ProperEulerAngles, AngleAxis]...}

# build expression for tuples of the parameterizations data
function expr_gen(typed, fields)
    expr = :(())  # build an expression for a tuple
    acc = typed ? typed_acc : untyped_acc
    append!(expr.args, [acc(field) for field in fields])
    return expr
end
untyped_acc(i) = isa(i, Int) ? :(X[$(i)]) : :(X.$(i))
typed_acc(i) = isa(i, Int) ? :(T(X[$(i)])) : :(T(X.$(i)))


#
# function to create a code block to convert to and from Vec and Vector rtypes
#
function add_methods(rot_type)

    #
    # Build expressions for building typed and untyped tuples using an Integer Index
    #
    untyped_indexed_expr = expr_gen(false, 1:numel(rot_type))
    typed_indexed_expr = expr_gen(true, 1:numel(rot_type))

    #
    # Build expressions for building typed and untyped tuples using a field
    #

    fields = fieldnames(rot_type)
    untyped_field_expr = (rot_type <: Mat) ? untyped_indexed_expr : expr_gen(false, fields[1:numel(rot_type)])
    typed_field_expr = (rot_type <: Mat) ? typed_indexed_expr : expr_gen(true, fields[1:numel(rot_type)])


    # the code block to return
    qb = quote end

    #
    # convert to / from a mutable vector
    #
    qn = quote

        # convert to a mutable vector
        @inline convert(::Type{Vector}, X::$(rot_type)) = vcat($(untyped_field_expr.args...))
        @inline convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(typed_field_expr.args...))
        @inline call{T <: Vector}(::Type{T}, X::$(rot_type)) = convert(T, X)

        # convert from a mutable vector
        @inline convert{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = T($(untyped_indexed_expr.args...))
        @inline call{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = convert(T, X)

    end
    append!(qb.args, qn.args)

    #
    # NaN checking
    #
    nan_expr = Expr(:||, :(isnan($(untyped_field_expr.args[1]))),  :(isnan($(untyped_field_expr.args[2]))))
    for i = 3:numel(rot_type)
        nan_expr = Expr(:||, nan_expr, :(isnan($(untyped_field_expr.args[i]))))
    end
    push!(qb.args, :( isnan(X::$(rot_type)) = $(nan_expr)))

    #
    # Sneak in an extra constructor that casts Ints to Floats (who wants Ints for these thing anyways)
    #
    int_cons1 = [:( $(symbol("x$(i)"))::Integer ) for i in 1:numel(rot_type)]
    int_cons2 = [:( Float64($(symbol("x$(i)"))) ) for i in 1:numel(rot_type)]
    int_expr = :(
                    call(::Type{$(rot_type)}, $(int_cons1...)) = call($(rot_type), $(int_cons2...))
                )
    push!(qb.args, int_expr)

    return qb
end

#
# Now go through and add everything
#
for rt in RotTypeList

    # println("type: $(rt), vars: $(numel(rt))")

    # add import / export
    eval(add_methods(rt))
end

