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

    # build expressions for the input and output tupples
    input_expr = :(())   # build an expression for a tuple
    output_expr = :(())  # build an expression for a tuple
    for i = 1:numel(rot_type)
        xi = symbol("x$(i)")
        push!(input_expr.args, :($(xi)::Int))
        push!(output_expr.args, :($(xi)))
    end

    
    # and build     
    if (n_params(rot_type) == 1)
        q = quote
        
            # a constructor in case all inputs where Ints (element types should be AbstractFloat)
            call(::Type{$(rot_type)}, $(input_expr.args...)) = $(rot_type){$(def_params[1])}($(output_expr.args...))

        end
    elseif (n_params(rot_type) == 2)
        q = quote

            # a constructor for the case where all inputs are Ints (element types should be AbstractFloat)
            call(::Type{$(rot_type)}, $(input_expr.args...)) = $(rot_type){$(def_params[1]), $(def_params[2])}($(output_expr.args...))
            call{T}(::Type{$(rot_type){T}}, $(input_expr.args...)) = $(rot_type){T, $(def_params[2])}($(output_expr.args...))
    
        end
    else
        q = quote end
    end
end


#
# function to create a code block to perform conversion to and from vector types
# 
function add_vector_conversions(rot_type)  # for 4 element representations

    # special case for building the mat type, we need a tuple for each column
    if (rot_type <: Mat)
        construct_expr = :(())
        for c in 1:size(rot_type,2)
            col_expr = :(())
            idx = (c-1)*size(rot_type,1) + (1:size(rot_type,1))
            append!(col_expr.args, [:(T(X[$(i)])) for i in idx])
            push!(construct_expr.args, col_expr)
        end

        # use indices not fieldnames to access members
        output_expr = :(())  # build an expression for a tuple
        append!(output_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])

    else
        # grab from input vectors using an index
        construct_expr = :(())
        append!(construct_expr.args, [:(T(X[$(i)])) for i in 1:numel(rot_type)])

        # use fieldnames to access members for exporting (I think its faster)
        fields = fieldnames(rot_type)
        output_expr = :(())  # build an expression for a tuple
        append!(output_expr.args, [:(T(X.$(fields[i]))) for i in 1:numel(rot_type)])

    end

    quote
    
        # convert to a fixed size a vector
        convert{T <: $(rot_type)}(::Type{Vec}, X::T) = convert(Vec{$(numel(rot_type)), eltype(T)}, X)
        convert{T <: Real}(::Type{Vec{$(numel(rot_type)),T}}, X::$(rot_type)) = Vec{$(numel(rot_type)),T}($(output_expr.args...))  
        if !($(rot_type) <: FixedArray)  # need this for quaternions
            call{T}(::Type{Vec{$(numel(rot_type)),T}}, X::$(rot_type)) = convert(Vec{$(numel(rot_type)),T}, X)
            call(::Type{Vec}, X::$(rot_type)) = convert(Vec, X)
        end

        # convert to mutable vector
        convert{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, X)
        convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat($(output_expr.args...))

        # convert from a fixed size a vector
        function convert{rT <: $(rot_type), U <: Real}(::Type{rT}, X::Vec{$(numel(rot_type)), U})
            oT = add_params(rT, U)                                                  # output type
            T = eltype(oT)                                                          # output element type
            oT($(construct_expr.args...))
        end
        
        # convert from mutable vector
        function convert{rT <: $(rot_type), U <: Real}(::Type{rT}, X::Vector{U})
            (length(X)) == $(numel(rot_type)) || error(@sprintf("Input vector should have length %i", $(numel(rot_type))))
            oT = add_params(rT, U)
            T = eltype(oT)
            oT($(construct_expr.args...))
        end

        # have to overload some call methods in FixedSizeArrays here...
        if ($(rot_type) <: FixedSizeArrays.FixedArray)  
            call{T <: $(rot_type), U <: Real}(::Type{T}, v::Vec{$(numel(rot_type)), U}) = convert(T, v)
            call{T <: $(rot_type), U <: Real}(::Type{T}, v::Vector{U}) = convert(T, v)      
        end

        # Quaternions has a weird constructor 3 element vector constructor 
        if ($(rot_type) == Quaternion)
            # luckily the below will be called in preference to call{T}(::Type{Quaternion}, X::Vector)        
            call{T}(::Type{Quaternion}, X::Vector{T}) = length(X) == 4 ? convert(Quaternion, X) : (length(X) == 3 ? Quaternion(0, X[1], X[2], X[3]) : error("Vector should have 3 or 4 elements"))
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
    append!(qb.args, Rotations.add_vector_conversions(rot_type).args)  

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



