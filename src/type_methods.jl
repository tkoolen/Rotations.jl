# defines common methods for the various rotation types


RotTypeList    = [RotMatrix, Quaternion, SpQuat, EulerAngles, ProperEulerAngles, AngleAxis]
Mat33Types     = [RotMatrix]
Vec4Types      = [Quaternion, AngleAxis]
Vec3Types      = [SpQuat, EulerAngles, ProperEulerAngles]


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
promote_type_sp(::Type{Any}, ::Type{Any}) = DefaultElType
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
        add_elltype(::Type{$(rot_type)}) = $(rot_type){$(def_element)}

        # do nothing when it is present
        add_elltype{T}(::Type{$(rot_type){T}}) = $(rot_type){T}

        # add the element type from another source when it's not present
        add_elltype{T <: Real}(::Type{$(rot_type)}, ::Type{T}) = $(rot_type){T}
        add_elltype{T}(::Type{$(rot_type)}, ::Type{T}) = $(rot_type){eltype(T)}

        # add the element type from another source when it is already present (ignore the other source)
        add_elltype{T,U}(::Type{$(rot_type){T}}, ::Type{U}) = $(rot_type){T}

        #
        # Stripping the element type
        #
        strip_elltype(::Type{$(rot_type)}) = $(rot_type)
        strip_elltype{T}(::Type{$(rot_type){T}}) = $(rot_type)

        #
        # Adding all parameters
        #
        add_params{T <: $(rot_type)}(::Type{T}) = add_elltype(T)                  # by itself
        add_params{T <: $(rot_type), U}(::Type{T}, ::Type{U}) = add_elltype(T, U) # from another source
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
        add_elltype{ORDER}(::Type{$(rot_type){ORDER}}) = $(rot_type){ORDER, $(def_element)}

        # do nothing when it is present
        add_elltype{ORDER, T}(::Type{$(rot_type){ORDER, T}}) = $(rot_type){ORDER, T}

        # add the element type from another source when it's not present 
        add_elltype{ORDER, T <: Real}(::Type{$(rot_type){ORDER}}, ::Type{T}) = $(rot_type){ORDER, T}
        add_elltype{ORDER, T}(::Type{$(rot_type){ORDER}}, ::Type{T}) = $(rot_type){ORDER, eltype(T)}

        # add the element type from another source when it is already present (ignore the other source)
        add_elltype{ORDER, T,U}(::Type{$(rot_type){ORDER, T}}, ::Type{U}) = $(rot_type){ORDER, T}

        #
        # Stripping the element type
        #
        strip_elltype(::Type{$(rot_type)}) = $(rot_type)
        strip_elltype{ORDER}(::Type{$(rot_type){ORDER}}) = $(rot_type){ORDER}
        strip_elltype{ORDER, T}(::Type{$(rot_type){ORDER, T}}) = $(rot_type){ORDER}

        #
        # Adding all parameters
        #
        add_params{T <: $(rot_type)}(::Type{T}) = add_elltype(add_order(T)) 
        add_params{T <: $(rot_type), U}(::Type{T}, ::Type{U}) = add_elltype(add_order(T), U)                   # eltype from another source
        add_params{T <: $(rot_type), U <: $(rot_type)}(::Type{T}, ::Type{U}) = add_elltype(add_order(T, U), U) # order and eltype from another source
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
# function to create a code block to perform conversion to and from vector types
# 
function conversions_V4(rot_type)  # for 4 element representations
    quote
    
        # convert to a fixed size a vector
        convert{T <: $(rot_type)}(::Type{Vec}, X::T) = convert(Vec{4, eltype(T)}, X)
        convert{T <: Real}(::Type{Vec{4,T}}, X::$(rot_type)) = Vec{4,T}(T(X[1]), T(X[2]), T(X[3]), T(X[4]))  
        if !($(rot_type) <: FixedArray)  # need this for quaternions
            call{T}(::Type{Vec{4,T}}, X::$(rot_type)) = convert(Vec{4,T}, X)
            call(::Type{Vec}, X::$(rot_type)) = convert(Vec, X)
        end

        # convert to mutable vector
        convert{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, X)
        convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat(T(X[1]), T(X[2]), T(X[3]), T(X[4]))

        # convert from a fixed size a vector
        function convert{T <: $(rot_type), U <: Real}(::Type{T}, v::Vec{4, U})
            oT = add_params(T,U)                        # output type
            elT = eltype(oT)                            # output element type
            oT(elT(v[1]), elT(v[2]), elT(v[3]), elT(v[4]))
        end
        
        # convert from mutable vector
        function convert{T <: $(rot_type), U <: Real}(::Type{T}, v::Vector{U})
            (length(v)) == 4 || error("Input vector should have length 4")
            oT = add_params(T, U)
            elT = eltype(oT)
            oT(elT(v[1]),elT(v[2]),elT(v[3]),elT(v[4]))
        end

        # Quaternions has a weird constructor 3 element vector constructor 
        if ($(rot_type) == Quaternion)
            # luckily the below will be called in preference to call{T}(::Type{Quaternion}, X::Vector)        
            call{T}(::Type{Quaternion}, X::Vector{T}) = length(X) == 4 ? convert(Quaternion, X) : (length(X) == 3 ? Quaternion(0, X[1], X[2], X[3]) : error("Vector should have 3 or 4 elements"))  
        end
    end
end

#
# function to crate a code block to perform conversion to and from vector types
# 
function conversions_V3(rot_type) # for 3 element representations
    quote

        # convert to a fixed size vector
        convert{T <: $(rot_type)}(::Type{Vec}, X::T) = convert(Vec{3, eltype(T)}, X)
        convert{T <: Real}(::Type{Vec{3,T}}, X::$(rot_type)) = Vec{3,T}(T(X[1]), T(X[2]), T(X[3]))
        

        # convert to mutable vector
        convert{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, X)
        convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = vcat(T(X[1]),T(X[2]),T(X[3]))

        # convert from a fixed size a vector
        function convert{T <: $(rot_type), U <: Real}(::Type{T}, v::Vec{3, U})
            oT = add_params(T, U)                       # output type
            elT = eltype(oT)                            # output element type
            oT(elT(v[1]), elT(v[2]), elT(v[3]))
        end
        
        # convert from mutable vector
        function convert{T <: $(rot_type), U <: Real}(::Type{T}, v::Vector{U})
            (length(v)) == 3 || error("Input vector should have length 3")
            oT = add_params(T, U)
            elT = eltype(oT)
            oT(elT(v[1]), elT(v[2]), elT(v[3]))
        end
    end
end


#
# function to crate a code block to perform conversion to and from vector types
# 
function conversions_Mat33(rot_type) # for 3 element representations
    quote

        # no need to add Vec conversion

        # convert to mutable vector
        convert{T <: $(rot_type)}(::Type{Vector}, X::T) = convert(Vector{eltype(T)}, vec(Matrix(X)))
        convert{T <: Real}(::Type{Vector{T}}, X::$(rot_type)) = convert(Vector{T}, vec(Matrix(X)))

        # convert from immutable vector
        function convert{T <: $(rot_type), U <: Real}(::Type{T}, v::Vec{9,U})
            oT = add_params(T, U)
            elT = eltype(oT)
            oT( (elT(v[1]), elT(v[2]), elT(v[3])), 
                (elT(v[4]), elT(v[5]), elT(v[6])), 
                (elT(v[7]), elT(v[8]), elT(v[9])) 
              )
        end
        call{T <: $(rot_type), U}(::Type{T}, v::Vec{9,U}) = convert(T, v)
       
        # convert from mutable vector
        function convert{T <: $(rot_type), U <: Real}(::Type{T}, v::Vector{U})
            (length(v)) == 9 || error("Input vector should have length 9")
            oT = add_params(T, U)
            elT = eltype(oT)
            convert(Mat{3,3, elT}, reshape(v, 3, 3))
        end
        call{T <: $(rot_type), U}(::Type{T}, v::Vector{U}) = convert(T, v)
    end
end

# the number of template parameters
n_params{T}(::Type{T}) = length(T.parameters)
n_params(::Type{RotMatrix}) = 1  # otherwise they'd be 3 because its a typecast ...


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
# add conversions to and from mutable and immutable vectors
# 
function add_vector_conversions(rot_type)
    if any(rot_type .== Vec4Types)
        conversions_V4(rot_type)
    elseif any(rot_type .== Vec3Types)
        conversions_V3(rot_type)
    elseif any(rot_type .== RotTypeList)
        conversions_Mat33(rot_type)
    end
end


#
# add NaN checking
# 
function add_nan_check(rot_type)
    if any(rot_type .== Vec4Types)
        quote isnan(X::$(rot_type))  = any(@fsa_isnan_vec(X, 4)) end
    elseif any(rot_type .== Vec3Types)
        quote isnan(X::$(rot_type))  = any(@fsa_isnan_vec(X, 3)) end
    elseif any(rot_type .== Mat33Types)
        quote isnan(X::$(rot_type))  = any(@fsa_isnan(X, 3, 3)) end
    end
end


#
# add a nvars function for the number of elements in the representations
# 
function add_nvars(rot_type)
    if any(rot_type .== Vec4Types)
        quote
            nvars(X::$(rot_type))  = nvars($(rot_type))
            nvars{T <: $(rot_type)}(::Type{T})  = 4
        end
    elseif any(rot_type .== Vec3Types)
        quote
            nvars(X::$(rot_type))  = nvars($(rot_type))
            nvars{T <: $(rot_type)}(::Type{T})  = 3
        end
    elseif any(rot_type .== Mat33Types)
        quote
            nvars(X::$(rot_type))  = nvars($(rot_type))
            nvars{T <: $(rot_type)}(::Type{T})  = 9
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

    # add the vector conversion code
    append!(qb.args, Rotations.add_vector_conversions(rot_type).args)  

    # the number of elements in the representation
    append!(qb.args, Rotations.add_nvars(rot_type).args)  

    # NaN checks
    append!(qb.args, Rotations.add_nan_check(rot_type).args)

    # return the code block
    return qb

end


#
# Now go through and add everything
# 
for t in RotTypeList
    #println(t)
    #println("n_params: $(n_params(t))")
    eval(add_methods(t))
end


# allow constuction from any rotation parameterization
call{T <: RotationTypes, U <: RotationTypes}(::Type{T}, X::U) = convert_rotation(add_params(T, U), X)



################################################################
# And some extra methods
################################################################

# add isanan to fixed size arrays (convenient)
@doc """
function to convert an immutable transformation matrix from RN to PN
"""  ->
function projective{T}(Tmat::RotMatrix{T})
    return @fsa_projective(Tmat, 3, 3, T)
end



