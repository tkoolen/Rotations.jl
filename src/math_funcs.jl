# define maths function for rotations types
import Base.(*), Base.(.*), Base.(-)

function add_maths(rot_type)

    # create an expression for a tuple X[1], X[2], etc...
    fields = fieldnames(rot_type)

    # the code block we are creating
    qb = quote; end

    # allow rotating a Vec{3, ...} using *    
    if !(rot_type <: FixedArray)
        qn = quote
            (*){T <: Real}(R::$(rot_type), X::Vec{3,T}) = rotate_point(R, X)
        end
        append!(qb.args, qn.args)
    end
    

    # allow multiplcation by a scalar (might make sense with small angle approximations)
    if !(rot_type <: FixedArray) && (rot_type != Quaternion)  

        # build the rhs epxression
        rhs_expr = :(())
        append!(rhs_expr.args, [:(s * X.$(field)) for field in fields[1:numel(rot_type)]])

        # form the whole expression
        qn = quote
            (*){T <: $(rot_type), U <: Real}(X::T, s::U) = strip_eltype(T)($(rhs_expr.args...))
            (*){T <: $(rot_type), U <: Real}(s::U, X::T) = strip_eltype(T)($(rhs_expr.args...))
        end
        append!(qb.args, qn.args)
    
    end


    # allow .* multiplcation by a Vec for thw same reason as me have multiplication by a scalar
    if !(rot_type <: FixedArray) && (rot_type != Quaternion)  

        # build the rhs epxression
        rhs_expr = :(())
        append!(rhs_expr.args, [:(s[$(i)] * X.$(fields[i])) for i in 1:numel(rot_type)])

        # form the whole expression
        qn = quote
            (.*){T <: $(rot_type), U <: Real}(X::T, s::Vec{$(numel(rot_type)), U}) = strip_eltype(T)($(rhs_expr.args...))
            (.*){T <: $(rot_type), U <: Real}(s::Vec{$(numel(rot_type)), U}, X::T) = strip_eltype(T)($(rhs_expr.args...))
        end
        append!(qb.args, qn.args)
    
    end

    # allow unary -
    if !(rot_type <: FixedArray) && (rot_type != Quaternion)  

        # build the rhs epxression
        rhs_expr = :(())
        append!(rhs_expr.args, [:(-X.$(fields[i])) for i in 1:numel(rot_type)])

        # form the whole expression
        qn = quote
            (-){T <: $(rot_type)}(X::T) = T($(rhs_expr.args...))
        end
        append!(qb.args, qn.args)
    
    end

    return qb
end







