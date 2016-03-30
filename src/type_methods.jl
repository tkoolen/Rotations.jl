# defines common methods for the various rotation types
RotTypeList    = [RotMatrix, Quaternion, SpQuat, EulerAngles, ProperEulerAngles, AngleAxis]

# It'd be nice if this was a abstract type...
RotationTypes  = Union{RotTypeList...}


#################################################################
# force the output type to have all of its template parameters
# N.B. these need to be generated for type stability because 
# T and U need to be known to the compiler to detemine the output type
#################################################################

# special promote methods ignoring "Any" when promoting types
promote_type_sp{T <: Real}(::Type{Any}, ::Type{T}) = T
promote_type_sp{T <: Real}(::Type{T}, ::Type{Any}) = T
promote_type_sp(::Type{Any}, ::Type{Any}) = @DefaultElType()
promote_type_sp{T <: Real, U <: Real}(::Type{T}, ::Type{U}) = promote_type(T, U)

# Ints aren't good for expressing rotations, but are handy to type
promote_type_sp{T <: Int}(::Type{Any}, ::Type{T}) = @DefaultElType()  # we need a float type
promote_type_sp{T <: Int}(::Type{T}, ::Type{Any}) = @DefaultElType()  # we need a float type
promote_type_sp{T <: Int, U <: Int}(::Type{T}, ::Type{U}) = @DefaultElType()
promote_type_sp{T <: Int, U <: Real}(::Type{T}, ::Type{U}) = U




#####################################################
# Generate the element type function before proceding
#####################################################

for rot_type in RotTypeList
    if n_params(rot_type) == 1
        qb = quote
            eltype(::Type{$(rot_type)}) = Any
            eltype{T}(::Type{$(rot_type){T}}) = T
            eltype{T <: $(rot_type)}(X::T) = eltype(T)
        end
    elseif n_params(rot_type) == 2
        qb = quote
            eltype(::Type{$(rot_type)}) = Any
            eltype{P1}(::Type{$(rot_type){P1}}) = Any
            eltype{P1, T}(::Type{$(rot_type){P1, T}}) = T
            eltype{T <: $(rot_type)}(X::T) = eltype(T)
        end
    end
    eval(qb)
end







#####################################################
# Worker macros
#####################################################

# and get index methods to each type
function add_indexing{rot_type}(::Type{rot_type})
    if rot_type <: FixedArray
        q = quote; end  # no need for subtypes of this
    else
        q = quote
            # access by index
            @inline getindex{T <: $(rot_type)}(x::T, i::Integer) = x.(i)
        end
    end
end

#
# function to crate a code block to add parameter stuff when the type has only the element type as a parameter
# 
function param_functions1{rot_type}(::Type{rot_type})
    def_element = default_params(rot_type)[end]
    quote
    
        #
        # Adding the element type
        #

        # add the element type when it's not present
        add_eltype(::Type{$(rot_type)}) = $(rot_type){$(def_element)}

        # do nothing when it is present
        add_eltype{T}(::Type{$(rot_type){T}}) = $(rot_type){T}

        # add the element type from another source when it's not present
        add_eltype{T <: AbstractFloat}(::Type{$(rot_type)}, ::Type{T}) = $(rot_type){T}
        add_eltype{T <: Real}(::Type{$(rot_type)}, ::Type{T}) = $(rot_type){$(def_element)}  # need a float
        add_eltype(::Type{$(rot_type)}, ::Type{Any}) = $(rot_type){$(def_element)}  # need a float
        add_eltype{T}(::Type{$(rot_type)}, ::Type{T}) = add_eltype($(rot_type), eltype(T))

        # add the element type from another source when it is already present (ignore the other source)
        add_eltype{T,U}(::Type{$(rot_type){T}}, ::Type{U}) = $(rot_type){T}

        #
        # Stripping the element type
        #
        strip_eltype(::Type{$(rot_type)}) = $(rot_type)
        strip_eltype{T}(::Type{$(rot_type){T}}) = $(rot_type)

        #
        # Adding all parameters
        #
        add_params{T <: $(rot_type)}(::Type{T}) = add_eltype(T)                  # by itself
        add_params{T <: $(rot_type), U}(::Type{T}, ::Type{U}) = add_eltype(T, U) # from another source
        add_params{T <: $(rot_type), U}(::Type{T}, X::U) = add_params(T, U)       # from another source

        # output type with a promoted element types (N.B. should I call this something else
        promote_eltype{T <: $(rot_type), U <: Real}(::Type{T}, ::Type{U}) = $(rot_type){promote_type_sp(eltype(T), U)}
        promote_eltype{T <: $(rot_type)}(::Type{T}, ::Type{Any}) = add_params(T)  
        function promote_eltype{T <: $(rot_type), U <: $(rot_type)}(::Type{T}, ::Type{U})
            oT = add_params(T, U)
            promote_eltype(oT, promote_type_sp(eltype(oT), eltype(U)))
        end
        promote_eltype{T <: $(rot_type), U}(::Type{T}, ::Type{U}) = $(rot_type){promote_type_sp(eltype(T), eltype(U))}

        # and do the actual promotion
        promote_eltype{T <: $(rot_type), U}(X::T, ::Type{U}) = convert(promote_eltype(T, U), X)  # via convert


    end
end


#
# function to crate a code block to add parameter stuff when the type has two parameters
# 
function param_functions2{rot_type}(::Type{rot_type})
    def_element = default_params(rot_type)[end]
    
    # get the default element type and order 
    def_order = default_params(rot_type)[1]
    def_element = default_params(rot_type)[end]

    order_type = super(def_order)
    quote
        
        #   
        # get the order
        #
        get_order(::Type{$(rot_type)}) = $(def_order)
        get_order{ORDER}(::Type{$(rot_type){ORDER}}) = ORDER
        get_order{ORDER,T}(::Type{$(rot_type){ORDER, T}}) = ORDER

        #
        # Adding the order
        #
        add_order(::Type{$(rot_type)}) = $(rot_type){$(def_order)}
        add_order{ORDER}(::Type{$(rot_type){ORDER}}) = $(rot_type){ORDER}
        add_order{ORDER,T}(::Type{$(rot_type){ORDER,T}}) = $(rot_type){ORDER,T}

        # add the order when it's not present from another source
        add_order{T <: $(rot_type)}(::Type{$(rot_type)}, ::Type{T}) = $(rot_type){get_order(T)}
        add_order{ORDER, T <: $(rot_type)}(::Type{$(rot_type){ORDER}}, ::Type{T}) = $(rot_type){ORDER}
        add_order{ORDER, T, U <: $(rot_type)}(::Type{$(rot_type){ORDER,T}}, ::Type{U}) = $(rot_type){ORDER,T}
        
        #
        # Dont strip the order
        #
    
        #
        # Adding the element type
        #

        # add the element type when it's not present
        add_eltype{ORDER}(::Type{$(rot_type){ORDER}}) = $(rot_type){ORDER, $(def_element)}

        # do nothing when it is present
        add_eltype{ORDER, T}(::Type{$(rot_type){ORDER, T}}) = $(rot_type){ORDER, T}

        # add the element type from another source when it's not present 
        add_eltype{ORDER, T <: AbstractFloat}(::Type{$(rot_type){ORDER}}, ::Type{T}) = $(rot_type){ORDER, T}
        add_eltype{ORDER, T <: Real}(::Type{$(rot_type){ORDER}}, ::Type{T}) = $(rot_type){ORDER, $(def_element)}  # need a float
        add_eltype{ORDER}(::Type{$(rot_type){ORDER}}, ::Type{Any}) = $(rot_type){ORDER, $(def_element)}  # need a float
        add_eltype{ORDER, T}(::Type{$(rot_type){ORDER}}, ::Type{T}) = add_eltype($(rot_type){ORDER}, eltype(T))

        # add the element type from another source when it is already present (ignore the other source)
        add_eltype{ORDER, T,U}(::Type{$(rot_type){ORDER, T}}, ::Type{U}) = $(rot_type){ORDER, T}

        #
        # Stripping the element type
        #
        strip_eltype(::Type{$(rot_type)}) = $(rot_type)
        strip_eltype{ORDER}(::Type{$(rot_type){ORDER}}) = $(rot_type){ORDER}
        strip_eltype{ORDER, T}(::Type{$(rot_type){ORDER, T}}) = $(rot_type){ORDER}

        #
        # Adding all parameters
        #
        add_params{T <: $(rot_type)}(::Type{T}) = add_eltype(add_order(T)) 
        add_params{T <: $(rot_type), U}(::Type{T}, ::Type{U}) = add_eltype(add_order(T), U)                   # eltype from another source
        add_params{T <: $(rot_type), U <: $(rot_type)}(::Type{T}, ::Type{U}) = add_eltype(add_order(T, U), U) # order and eltype from another source
        add_params{T <: $(rot_type), U}(::Type{T}, X::U) = add_params(T, U)       # from another source

        # output type with a promoted element types
        promote_eltype{T <: $(rot_type), U <: Real}(::Type{T}, ::Type{U}) = $(rot_type){get_order(T), promote_type_sp(eltype(T), U)}
        promote_eltype{T <: $(rot_type)}(::Type{T}, ::Type{Any}) = add_params(T)  
        function promote_eltype{T <: $(rot_type), U <: $(rot_type)}(::Type{T}, ::Type{U})
            oT = add_params(T, U)
            promote_eltype(oT, promote_type_sp(eltype(oT), eltype(U)))
        end
        promote_eltype{T <: $(rot_type), U}(::Type{T}, ::Type{U}) = promote_eltype(add_params(T), eltype(U))


        # and do the actual promotion
        promote_eltype{T <: $(rot_type), U}(X::T, ::Type{U}) = convert(promote_eltype(T, U), X)  

    end
end

#
# function to add constructors for each type
# 
function add_constructors(rot_type)

    def_params = default_params(rot_type)  # default element to use

    # build expressions for the input / output elements
    xsym = [symbol("x$(i)") for i in 1:numel(rot_type)]

    # create an expression for a tuple x1, x2, x3...
    rhs_expr = :(())
    append!(rhs_expr.args, [:($(xsym[i])) for i in 1:numel(rot_type)])

    # create an expression for a tuple X[1], X[2], etc...
    if !(rot_type <: Mat)
        rhs_tuple_expr = :(())
        append!(rhs_tuple_expr.args, [:(X[$(i)]) for i in 1:numel(rot_type)])

        # for element conversion
        fields = fieldnames(rot_type)
        rhs_typed_fields_expr = :(())
        append!(rhs_typed_fields_expr.args, [:(T(X.$(fields[i]))) for i in 1:numel(rot_type)])

    else

        # I'm not sure what's going on with the matrix constructor, but its not what I want        
        rhs_tuple_expr = :(())
        for c in 1:size(rot_type,2)
            col_expr = :(())
            idx = (c-1)*size(rot_type,1) + (1:size(rot_type,1))
            append!(col_expr.args, [:(X[$(i)]) for i in idx])
            push!(rhs_tuple_expr.args, col_expr)
        end
    end

    # create an expression for a tuple x1::Int, x2::Int, x3::Int... for allowing construction from Ints
    lhs_expr = :(())
    append!(lhs_expr.args, [:($(xsym[i])::T) for i in 1:numel(rot_type)])

    #
    # add extra construction methods for constructing from NTupples and FixedSizeArrays
    #
    qb = quote

        # allow construction from a tuple if this isn't a fixed array
        convert{T <: $(rot_type), U <: Real}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = add_params(T, U)($(rhs_tuple_expr.args...))

        # allow construction from an Fixed Size Array Vector
        convert{T <: $(rot_type), U <: Real}(::Type{T}, X::Vec{$(numel(rot_type)), U}) = add_params(T, U)($(rhs_tuple_expr.args...))

    end
    
    # realign call and convert if this was a fixed size array
    if (rot_type <: Mat)
        qn = quote
            call{T <: $(rot_type), U <: Real}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = convert(T, X)
            call{T <: $(rot_type), U <: Real}(::Type{T}, X::Vec{$(numel(rot_type)), U}) = convert(T, X)

        end
        append!(qb.args, qn.args)
    end

    #
    # add code for creating from a mutable Vector
    #
    if (rot_type == Quaternion)
        qn = quote

            # Quaternions come with a Vector conversion defined, but not the one we want
            function convert{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) 
                (length(X) == 4) ? add_params(T,U)($(rhs_tuple_expr.args...)) : (length(X) == 3) ? add_params(T,U)(0, X[1], X[2], X[3]) : error("Vector should have 3 or 4 elements")
            end
            call{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = convert(T, X)

        end
    elseif (rot_type <: Mat)
        qn = quote
            convert{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = add_params(T, U)($(rhs_tuple_expr.args...))
            call{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = convert(T, X)
        end
    else
        qn = quote
            convert{T <: $(rot_type), U <: Real}(::Type{T}, X::Vector{U}) = add_params(T, U)($(rhs_tuple_expr.args...))  
        end
    end
    append!(qb.args, qn.args)

    
    #
    # add element conversion methods
    # N.B. FSA's and Quaternions already have element conversion
    #
    if (rot_type != Quaternion) && (!(rot_type <: Mat))
        qn = quote
            convert(::Type{$(rot_type)}, X::$(rot_type)) = X
            convert{T <: $(rot_type)}(::Type{T}, X::T) = X
            if ($(n_params(rot_type)) == 1)
                convert{T <: Real, U <: Real}(::Type{$(rot_type){T}}, X::$(rot_type){U}) = $(rot_type){T}($(rhs_typed_fields_expr.args...))
            end
        end
        append!(qb.args, qn.args)
    end
    

    #
    # add extra construction methods for constructing from Integer Vectors and FixedSizeArrays (we dont want ints for rotations)
    #
    if !(rot_type <: Mat)
        qn = quote

            # allow construction from integer inputs 
            call{T <: Integer}(::Type{$(rot_type)}, $(lhs_expr.args...)) = add_params($(rot_type))($(rhs_expr.args...))

            # allow construction from an integer tuple
            convert{rT <: $(rot_type), T <: Integer}(::Type{rT}, X::NTuple{$(numel(rot_type)), T}) = add_params(rT)($(rhs_tuple_expr.args...))            

            # allow construction from an Fixed Size Array Vector of ints
            convert{rT <: $(rot_type), T <: Integer}(::Type{rT}, X::Vec{$(numel(rot_type)), T}) = add_params(rT)($(rhs_tuple_expr.args...))

            # Quaternions has a weird constructor 3 element vector constructor 
            if ($(rot_type) == Quaternion)
                function convert{rT <: $(rot_type), T <: Integer}(::Type{rT}, X::Vector{T}) 
                    (length(X) == 4) ? add_params($(rot_type))($(rhs_tuple_expr.args...)) : (length(X) == 3) ? add_params($(rot_type))(0, X[1], X[2], X[3]) : error("Vector should have 3 or 4 elements")
                end
            else
                # allow construction from an mutable Vector of Ints
                function convert{rT <: $(rot_type), T <: Integer}(::Type{rT}, X::Vector{T})
                    length(X) == $(numel(rot_type)) || error("Vector should have $(numel(rot_type)) elements")
                    add_params(rT)($(rhs_tuple_expr.args...))
                end
            end

        end
        append!(qb.args, qn.args)
    end

    # an extra methods if the type has two parameters
    if (n_params(rot_type) == 2)
        P1Type = super(def_params[1])
        qn = quote
            
            # and extra constructors
            call{T <: Real}(::Type{$(rot_type)}, $(lhs_expr.args...)) = add_params($(rot_type), T)($(rhs_expr.args...))
            call{P1 <: $(P1Type), T <: Real}(::Type{$(rot_type){P1}}, $(lhs_expr.args...)) = $(rot_type){P1, T}($(rhs_expr.args...))

            call{T <: Int}(::Type{$(rot_type)}, $(lhs_expr.args...)) = add_params($(rot_type))($(rhs_expr.args...))
            call{P1 <: $(P1Type), T <: Int}(::Type{$(rot_type){P1}}, $(lhs_expr.args...)) = add_params($(rot_type){P1})($(rhs_expr.args...))
    
            # extra element conversion
            convert{P1 <: $(P1Type), T <: Real}(::Type{$(rot_type){P1}}, X::$(rot_type){P1, T}) = X
            convert{P1 <: $(P1Type), T <: Real, U <: Real}(::Type{$(rot_type){P1,T}}, X::$(rot_type){P1,U}) = $(rot_type){P1, T}($(rhs_typed_fields_expr.args...))
        end
        append!(qb.args, qn.args)
    end

    return qb
    
end


#
# function to create a code block to convert to Vec and Vectro types
# 
function add_export_conversions(rot_type)  

    # special case for building the mat type, we need a tuple for each column
    if (rot_type <: Mat)
        # use indices not fieldnames to access members
        
        # untyped version
        output_expr = :(())  # build an expression for a tuple
        append!(output_expr.args, [:(X[$(i)]) for i in 1:numel(rot_type)])

        # typed version
        typed_output_expr = :(())  # build an expression for a tuple
        append!(typed_output_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])
    else

        # use fieldnames to access members for exporting (I think its faster)
        fields = fieldnames(rot_type)

        # untyped version
        output_expr = :(())  # build an expression for a tuple
        append!(output_expr.args, [:(X.$(fields[i])) for i in 1:numel(rot_type)])
        
        # typed version
        typed_output_expr = :(())  # build an expression for a tuple
        append!(typed_output_expr.args, [:(T(X.$(fields[i]))) for i in 1:numel(rot_type)])

    end

    # grab from input vectors using an index
    construct_expr = :(())
    append!(construct_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])

    # the code block to return
    qb = quote end

    if !(rot_type <: FixedArray)
        qn = quote

            # convert to a fixed size a vector
            convert{T <: $(rot_type)}(::Type{Vec}, X::T) = Vec{$(numel(rot_type)), eltype(T)}($(output_expr.args...))
            convert{T <: Real, rT <: $(rot_type)}(::Type{Vec{$(numel(rot_type)), T}}, X::rT) = Vec{$(numel(rot_type)), T}($(typed_output_expr.args...))
            call{T <: Vec, rT <: $(rot_type)}(::Type{T}, X::rT) = convert(T, X)

            # convert to mutable vector
            convert{T <: $(rot_type)}(::Type{Vector}, X::T) = vcat($(output_expr.args...))
            convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(typed_output_expr.args...))

        end
        append!(qb.args, qn.args)
    elseif (rot_type <: Mat)
        qn = quote
            # convert to mutable vector
            convert{T <: $(rot_type)}(::Type{Vector}, X::T) = vcat($(output_expr.args...))                # appears to be broken in FixedSizeArrays v0.0.9
            convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(typed_output_expr.args...)) 
            call{T <: Vector}(::Type{T}, X::$(rot_type)) = convert(T, X)
        end
        append!(qb.args, qn.args)
    end
    

    # convert to an ntuple
    qn = quote
        convert{T <: $(rot_type)}(::Type{Tuple}, X::T) = $(output_expr)
        convert{T <: $(rot_type)}(::Type{NTuple}, X::T) = convert(NTuple{$(numel(rot_type)), eltype(T)}, X)
        convert{T <: Real, rT <: $(rot_type)}(::Type{NTuple{$(numel(rot_type)), T}}, X::rT) = $(typed_output_expr)
    end
    append!(qb.args, qn.args)
    return qb
end



#
# add parameter filling / copying etc
# 
function add_param_functions(rot_type)
    if (Rotations.n_params(rot_type) == 1)
        Rotations.param_functions1(rot_type)
    elseif (Rotations.n_params(rot_type) == 2)
        Rotations.param_functions2(rot_type)
    else
        error("param_functions: macro not defined for types with more than two parameters")
    end
end


#
# add NaN checking
# 
function add_nan_check(rot_type)

    if (rot_type <: Mat)
        # build a big or expression
        output_expr = Expr(:||, :(isnan(X[$(numel(rot_type) - 1)])), :(isnan(X[$(numel(rot_type))])))
        for i in numel(rot_type)-2:-1:1
            output_expr = Expr(:||,  :(isnan(X[$(i)])), output_expr)
        end
    else

        # use fieldnames to access members for exporting (I think its faster)
        fields = fieldnames(rot_type)

        # build a big or expression
        output_expr = Expr(:||, :(isnan(X.$(fields[numel(rot_type) - 1]))), :(isnan(X.$(fields[numel(rot_type)]))))
        for i in numel(rot_type)-2:-1:1
            output_expr = Expr(:||,  :(isnan(X.$(fields[i]))) , output_expr)
        end
    end
    quote 
        isnan(X::$(rot_type))  = $(output_expr)
    end
end



#
# Add all methods for a specific type
# 
include("math_funcs.jl")  # add some maths functions as well
function add_methods(rot_type)

    qb = quote end

    # add indexing scheme
    append!(qb.args, Rotations.add_indexing(rot_type).args)

    # add the parameter manipulation methods
    append!(qb.args, Rotations.add_param_functions(rot_type).args)

    # add extra constructors
    append!(qb.args, Rotations.add_constructors(rot_type).args)    

    # add the vector conversion code
    append!(qb.args, Rotations.add_export_conversions(rot_type).args)  

    # add maths
    append!(qb.args, Rotations.add_maths(rot_type).args)

    # NaN checks
    append!(qb.args, Rotations.add_nan_check(rot_type).args)

    # return the code block
    return qb

end


#
# Now go through and add everything
# 
for rt in RotTypeList
    #println(rt)
    #println("type: $(rt), params: $(n_params(rt)), vars: $(numel(rt))")
    eval(add_methods(rt))
end


# allow constuction from any rotation parameterization
call{T <: RotationTypes, U <: RotationTypes}(::Type{T}, X::U) = convert_rotation(add_params(T, U), X)



################################################################
# And some extra methods
################################################################

#@doc """
#function to convert an immutable transformation matrix from RN to PN
#"""  ->
#function projective{T}(Tmat::RotMatrix{T})
#    return @fsa_projective(Tmat, 3, 3, T)
#end



