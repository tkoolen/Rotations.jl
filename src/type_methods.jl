# defines common methods for the various rotation types
RotTypeList    = [RotMatrix, Quaternion, SpQuat, EulerAngles, ProperEulerAngles, AngleAxis]

# It'd be nice if this was a abstract type...
RotationTypes  = Union{RotTypeList...}

# function to get the number of elements in each parameterization
numel{T <: RotationTypes}(::Type{T}) = length(fieldnames(T))
numel(::Type{RotMatrix}) = 9 # special case
numel(::Type{Quaternion}) = 4 # Quaternions have an extra bool field

# the number of template parameters
n_params{T <: RotationTypes}(::Type{T}) = length(T.parameters)
n_params(::Type{RotMatrix})  = 1 # The behaviour of length(RotMatrix.parameters) seems unstable...



#################################################################
# force the output type to have all of its template parameters
# N.B. these need to be generated for type stability because 
# T and U need to be known to the compiler to detemine the output type
#################################################################

# special promote methods ignoring "Any" when promoting types
promote_type_sp{T <: AbstractFloat}(::Type{Any}, ::Type{T}) = T
promote_type_sp{T <: Real}(::Type{Any}, ::Type{T}) = @DefaultElType()  # we need a float type

promote_type_sp{T <: AbstractFloat}(::Type{T}, ::Type{Any}) = T
promote_type_sp{T <: Real}(::Type{T}, ::Type{Any}) = @DefaultElType()  # we need a float type

promote_type_sp(::Type{Any}, ::Type{Any}) = @DefaultElType()
promote_type_sp{T <: Real, U <: Real}(::Type{T}, ::Type{U}) = promote_type(T, U)



#####################################################
# Worker macros
#####################################################

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
    
    # and build     
    q = quote

        #
        # Allow construction from Ints by converting them
        #
    
        # allow construction from integer inputs
        if ($(n_params(rot_type)) == 1)  # need a special version for mor than 1 parameter
            call{T <: Integer}(::Type{$(rot_type)}, $(lhs_expr.args...)) = add_params($(rot_type))($(rhs_expr.args...))
        end

        # allow construction from an integer tuple
        convert{T <: Integer}(::Type{$(rot_type)}, X::NTuple{$(numel(rot_type)), T}) = add_params($(rot_type))($(rhs_tuple_expr.args...))            
        call{T <: Integer}(::Type{$(rot_type)}, X::NTuple{$(numel(rot_type)), T}) = add_params($(rot_type))($(rhs_tuple_expr.args...))

        # allow construction from an Fixed Size Array Vector of ints
        convert{T <: Integer}(::Type{$(rot_type)}, X::Vec{$(numel(rot_type)), T}) = add_params($(rot_type))($(rhs_tuple_expr.args...))
        call{T <: Integer}(::Type{$(rot_type)}, X::Vec{$(numel(rot_type)), T}) = add_params($(rot_type))($(rhs_tuple_expr.args...))

        # Quaternions has a weird constructor 3 element vector constructor 
        if ($(rot_type) == Quaternion)
            function convert{T <: Integer}(::Type{Quaternion}, X::Vector{T}) 
                (length(X) == 4) ? add_params($(rot_type))($(rhs_tuple_expr.args...)) : (length(X) == 3) ? add_params($(rot_type))(0, X[1], X[2], X[3]) : error("Vector should have 3 or 4 elements")
            end
        else
            # allow construction from an mutable Vector of Ints
            function convert{T <: Integer}(::Type{$(rot_type)}, X::Vector{T})
                length(X) == $(numel(rot_type)) || error("Vector should have $(numel(rot_type)) elements")
                add_params($(rot_type))($(rhs_tuple_expr.args...))
            end
        end

        #
        # Allow extra constructors to non-fixed array types (and fix the rotation matrix ones zzz)
        #
        if !($(rot_type) <: FixedArray) || ($(rot_type) <: Mat)

            # allow construction from a tuple if this isn't a fixed array
            convert{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = add_params(T, U)($(rhs_tuple_expr.args...))
            call{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = convert(add_params(T, U), X)

            # allow construction from an Fixed Size Array Vector
            convert{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vec{$(numel(rot_type)), U}) = add_params(T, U)($(rhs_tuple_expr.args...))
            call{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vec{$(numel(rot_type)), U}) = convert(add_params(T, U), X)

            # special case for Quaternion from a Vector
            if ($(rot_type) == Quaternion)
                function convert{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vector{U}) 
                    (length(X) == 4) ? add_params(T,U)($(rhs_tuple_expr.args...)) : (length(X) == 3) ? add_params(T,U)(0, X[1], X[2], X[3]) : error("Vector should have 3 or 4 elements")
                end
                call{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vector{U}) = convert(T, X)
            else
                convert{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vector{U}) = add_params(T, U)($(rhs_tuple_expr.args...))
                call{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vector{U}) = convert(add_params(T, U), X)
            end
        end
    end

    #
    # add the extra parameter methods if it has an extra parameter
    #

    if (n_params(rot_type) == 2)

        # create an expression for a tuple T(x1), T(x2), T(x3)...
        rhs_typed_expr = :(())
        append!(rhs_typed_expr.args, [:(T($(xsym[i]))) for i in 1:numel(rot_type)])

        rhs_tuple_typed_expr = :(())
        append!(rhs_tuple_typed_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])

        # the type for the first parameter
        P1Type = super(def_params[1])

        qn = quote

            #
            # N.B. with the extra parameter, the fixed size array constructors seems to need spoon feeding...
            #
            
            call{T <: AbstractFloat}(::Type{$(rot_type)}, $(lhs_expr.args...)) = add_params($(rot_type), T)($(rhs_expr.args...))
            call{P1 <: $(P1Type), T <: AbstractFloat}(::Type{$(rot_type){P1}}, $(lhs_expr.args...)) = add_params($(rot_type){P1}, T)($(rhs_expr.args...))

            # define element conversion
            convert{P1 <: $(P1Type), T <: AbstractFloat}(::Type{$(rot_type)}, X::$(rot_type){P1, T}) = X
            convert{P1 <: $(P1Type), T <: AbstractFloat}(::Type{$(rot_type){P1}}, X::$(rot_type){P1, T}) = X
            convert{P1 <: $(P1Type), T <: AbstractFloat}(::Type{$(rot_type){P1, T}}, X::$(rot_type){P1, T}) = X
            convert{P1 <: $(P1Type), T <: AbstractFloat, U <: AbstractFloat}(::Type{$(rot_type){P1, T}}, X::$(rot_type){P1, U}) = $(rot_type){P1, T}($(rhs_tuple_typed_expr.args...))
            
            call{T <: AbstractFloat}(::Type{$(rot_type)}, X::NTuple{$(numel(rot_type)), T}) = add_params($(rot_type), T)($(rhs_tuple_expr.args...))
            call{P1 <: $(P1Type), T <: AbstractFloat}(::Type{$(rot_type){P1}}, X::NTuple{$(numel(rot_type)), T}) = add_params($(rot_type){P1}, T)($(rhs_tuple_expr.args...))
            convert{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::NTuple{$(numel(rot_type)), U}) = add_params(T, U)($(rhs_tuple_expr.args...))
            
            call{T <: AbstractFloat}(::Type{$(rot_type)}, X::Vec{$(numel(rot_type)), T}) = add_params($(rot_type), T)($(rhs_tuple_expr.args...))
            call{P1 <: $(P1Type), T <: AbstractFloat}(::Type{$(rot_type){P1}}, X::Vec{$(numel(rot_type)), T}) = add_params($(rot_type){P1}, T)($(rhs_tuple_expr.args...))
            convert{T <: $(rot_type), U <: AbstractFloat}(::Type{T}, X::Vec{$(numel(rot_type)), U}) = add_params(T, U)($(rhs_tuple_expr.args...))            

            # allow construction from integer inputs
            call{P1, T <: Integer}(::Type{$(rot_type){P1}}, $(lhs_expr.args...)) = add_params($(rot_type){P1})($(rhs_expr.args...))
            
            convert{P1, T <: Integer}(::Type{$(rot_type){P1}}, X::NTuple{$(numel(rot_type)), T}) = add_params($(rot_type){P1})($(rhs_tuple_expr.args...))
            call{P1, T <: Integer}(::Type{$(rot_type){P1}}, X::NTuple{$(numel(rot_type)), T}) = convert($(rot_type){P1}, X) 

            convert{P1, T <: Integer}(::Type{$(rot_type){P1}}, X::Vec{$(numel(rot_type)), T}) = add_params($(rot_type){P1})($(rhs_tuple_expr.args...))    
            call{P1, T <: Integer}(::Type{$(rot_type){P1}}, X::Vec{$(numel(rot_type)), T}) = convert($(rot_type){P1}, X)

        end
        append!(q.args, qn.args)
    end

    

    return q
    
end


#
# function to create a code block to convert to Vec and Vectro types
# 
function add_export_conversions(rot_type)  

    # special case for building the mat type, we need a tuple for each column
    if (rot_type <: Mat)
        # use indices not fieldnames to access members
        output_expr = :(())  # build an expression for a tuple
        append!(output_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])
    else

        # use fieldnames to access members for exporting (I think its faster)
        fields = fieldnames(rot_type)
        output_expr = :(())  # build an expression for a tuple
        append!(output_expr.args, [:(T(X.$(fields[i]))) for i in 1:numel(rot_type)])

    end

    # grab from input vectors using an index
    construct_expr = :(())
    append!(construct_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])

    quote

        if !($(rot_type) <: FixedArray)  # fixed arrays have these exports defined
    
            # convert to a fixed size a vector
            convert{T <: $(rot_type)}(::Type{Vec}, X::T) = convert(Vec{$(numel(rot_type)), eltype(T)}, X)
            convert{T <: Real}(::Type{Vec{$(numel(rot_type)),T}}, X::$(rot_type)) = Vec{$(numel(rot_type)),T}($(output_expr.args...))  
            call{T <: Vec, U <:  $(rot_type)}(::Type{T}, X::U) = convert(T, X)
            

            # convert to mutable vector
            convert{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, X)
            convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(output_expr.args...))

        else
            call{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(output_expr.args...)) 
        end

        

        # convert to an ntuple
        convert{T <: $(rot_type)}(::Type{Tuple}, X::T) = convert(NTuple{$(numel(rot_type)), eltype(T)}, X)
        convert{T <: $(rot_type)}(::Type{NTuple}, X::T) = convert(NTuple{$(numel(rot_type)), eltype(T)}, X)
        convert{T <: Real}(::Type{NTuple{$(numel(rot_type)), T}}, X::T) = $(output_expr.args)

        # This can't be good...
        if ($(rot_type) <: Mat)

            # convert to a fixed size vector    
            convert{T <: $(rot_type)}(X::T, ::Type{Vec}) = convert(Vec{$(numel(rot_type)), eltype(T)}, X)
            convert{T <: Real}(::Type{Vec{$(numel(rot_type)),T}}, X::$(rot_type)) = Vec{$(numel(rot_type)),T}($(output_expr.args...))
            call{T <: $(rot_type)}(::Type{Vec}, X::T) = convert(Vec{$(numel(rot_type)), eltype(T)}, X)
            call{T <: Real}(::Type{Vec{$(numel(rot_type)),T}}, X::$(rot_type)) = convert(Vec{$(numel(rot_type)),T}, X)

            # convert to a mutable vector    
            convert{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, X)
            convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(output_expr.args...))
            call{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, X)
            call{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = convert(Vector{T}, X)

        end

    end
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
    if (rot_type <: FixedSizeArrays.Mat)  # special case
        quote 
            isnan(X::$(rot_type))  = any(@fsa_isnan(X, $(size(rot_type,1)), $(size(rot_type,2)))) 
        end
    else
        quote 
            isnan(X::$(rot_type))  = any(@fsa_isnan_vec(X, $(numel(rot_type)))) 
        end
    end
end



#
# Add all methods for a specific type
# 
function add_methods(rot_type)

    qb = quote end

    # add the parameter manipulation methods
    append!(qb.args, Rotations.add_param_functions(rot_type).args)

    # add the parameter manipulation methods
    append!(qb.args, Rotations.add_constructors(rot_type).args)    

    # add the vector conversion code
    append!(qb.args, Rotations.add_export_conversions(rot_type).args)  

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

@doc """
function to convert an immutable transformation matrix from RN to PN
"""  ->
function projective{T}(Tmat::RotMatrix{T})
    return @fsa_projective(Tmat, 3, 3, T)
end



