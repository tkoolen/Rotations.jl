# function to perform tests of the rotation functions in the Rotations module
using Base.Test
using FixedSizeArrays
using Quaternions
using Rotations

verbose_rot_test = false


numel{T <: Rotations.RotationTypes}(X::T) = Rotations.numel(Rotations.strip_eltype(T))
numel{T}(X::T) = (T <: FixedSizeArrays.FixedArray) || (T <: FixedSizeArrays.AbstractArray)  ? length(X) : error("numel undefined for input of type $(T)") 

# a macro to test if the contents of two containers are approximatley equal
macro contents_approx_eq(a, b)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test typeof(ai) == typeof(bi)
            @test_approx_eq ai bi
        end
    end
end
macro contents_approx_eq_eps(a, b, eps)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test typeof(ai) == typeof(bi)
            @test_approx_eq_eps ai bi $(esc(eps))
        end
    end
end


# a macro to test if tow types are approximatey equal
macro types_approx_eq(a, b)
    quote
        @test typeof($(esc(a))) == typeof($(esc(b)))
        @contents_approx_eq($(esc(a)), $(esc(b)))
    end
end
macro types_approx_eq_eps(a, b, eps)
    quote
        @test typeof($(esc(a))) == typeof($(esc(b)))
        @contents_approx_eq_eps($(esc(a)), $(esc(b)), $(esc(eps)))
    end
end

# a macro to test if the contents of two containers are approximatley equal, without examining the element types
macro contents_approx_eq_notype(a, b)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test_approx_eq ai bi
        end
    end
end


# build a full list of rotation types including the different ordering schemas
rot_types = Vector{Any}(0)
for rT in Rotations.RotTypeList
    if (Rotations.n_params(rT) == 2)

        # get the super type for the order parameter
        order_type = super(Rotations.default_params(rT)[1])
        for order in subtypes(order_type)
            push!(rot_types, rT{order})
        end
    end

    # ordered ones should have defaults so leave them in
    push!(rot_types, rT)
end

# define null rotations for conveniences
null_rotation(::Type{RotMatrix}) = RotMatrix(eye(3))
null_rotation(::Type{Quaternion}) = Quaternion(1.0,0,0,0)
null_rotation(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0, 0)
null_rotation(::Type{SpQuat}) = SpQuat(0.0,0.0,0.0)
null_rotation(::Type{EulerAngles}) = EulerAngles(0.0,0.0,0.0)
null_rotation(::Type{ProperEulerAngles}) = ProperEulerAngles(0.0,0.0,0.0)
null_rotation{ORD}(::Type{EulerAngles{ORD}}) = EulerAngles{ORD, Float64}(0.0,0.0,0.0)
null_rotation{ORD}(::Type{ProperEulerAngles{ORD}}) = ProperEulerAngles{ORD, Float64}(0.0,0.0,0.0)


#########################################################################
# Check fixed relationships
#########################################################################

# Do no rotation
R = RotMatrix(eye(3))
verbose_rot_test ? println("********************************\nIndentity checks\n********************************\n") : nothing
for rT in rot_types
    verbose_rot_test ? println(rT) : nothing
    rot_var = rT(R)
    null_var = null_rotation(rT)
    @types_approx_eq(rot_var, null_var)
end

#########################################################################
# Check conversion to and from mutable and immutable vectors,
# as well as element conversions
#########################################################################


verbose_rot_test ? println("\n\n\n********************************\nVector conversion checks\n********************************\n") :  nothing
R = RotMatrix(eye(3))
eltypes = subtypes(AbstractFloat)  # test different AbstractFloats
for rT in rot_types

    verbose_rot_test ? println("$(rT)") : nothing
    rot_var = rT(R)
    
    # export to immutable
    ivu = Vec(rot_var)
    @contents_approx_eq(rot_var, ivu)

    # export to mutable
    mvu = Vector(rot_var)
    @contents_approx_eq(rot_var, mvu)

    # import from immutable
    rot_ivu = rT(ivu)
    @types_approx_eq(rot_var, rot_ivu)

    # import from mutable
    rot_ivu = rT(mvu)
    @types_approx_eq(rot_var, rot_ivu)

    # test typed stuff
    for eT in eltypes
        
        #println("$(rT): $(eT)")
        
        # export to immutable
        ivt = Vec{Rotations.numel(rT), eT}(rot_var)
        ivc = convert(Vec{Rotations.numel(rT), eT}, ivu)
        @types_approx_eq(ivt, ivc)
    
        # export to mutable
        mvt = Vector{eT}(rot_var)
        mvc = convert(Vector{eT}, mvu)
        @types_approx_eq(mvt, mvc)

        # import from immutable
        rot_ivt = rT(ivt)
        @contents_approx_eq_notype(rot_var, rot_ivt)

        # import from mutable
        rot_mvt = rT(mvt)
        @contents_approx_eq_notype(rot_var, rot_mvt)

        # and test the element conversion on the rotation parameterization directly
        rot_c = convert(Rotations.add_params(rT, eT), rot_var)
        @contents_approx_eq(rot_c, rot_mvt)

    end
end


#########################################################################
# Test conversions between rotation types
#########################################################################

verbose_rot_test ? println("\n\n\n********************************\nRotation conversion checks\n********************************\n") : nothing

# test random round trip conversions
repeats = 1000
thresh = 1e-6
eye3 = @fsa([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
srand(0)
for rT_in in rot_types

    for rT_out in rot_types

        verbose_rot_test ? println("$(rT_in) - > $(rT_out)") : nothing

        # and each test
        #fcount = 0
        for i = 1:repeats
            
            # start with a random quaternion
            q = nquatrand()
            X = convert_rotation(rT_in, q)

            # round trip conversion
            Xd = convert_rotation(rT_in, convert_rotation(rT_out, X))

            # compare rotations before and after the round trip
            Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the identity
            rd = vecnorm(eye3 - Rout)
            @test rd <= thresh
            #fcount += (rd > thresh)
        end
        #if (fcount > 0)
        #    warn("Failed $(fcount) / $(repeats): $(rT_in) - > $(rT_out)")
        #end
    end
end



#########################################################################
# Actually rotate some stuff
#########################################################################

verbose_rot_test ? println("\n\n\n********************************\nTesting Rotations\n********************************\n") : nothing

# a random rotation
repeats = 1000
srand(0)

# rotate using each different parameterization
for rT in rot_types
    
    verbose_rot_test ? println("Rotating using the $(rT) parameterization") : nothing

    for i = 1:repeats

        Rm = RotMatrix(nquatrand())
        X = Vec{3, Float64}(randn(), randn(), randn())  # a point to rotate
        Xo = Rm * X  # Fixed size arrays better get this right
    
        R = rT(Rm)  # convert R to this formulation
        Xo_t = R * X
        @types_approx_eq_eps(Xo_t, Xo, 1e-10)

    end
end



























