# function to perform tests of the rotation functions in the Rotations module
using Base.Test
using FixedSizeArrays
using Quaternions
using Rotations

verbose_rot_test = false

import Rotations: numel
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
    if (length(rT.parameters) == 2)

        # get the super type for the order parameter
        order_type = super(Rotations.euler_order(rT))
        for order in subtypes(order_type)
            push!(rot_types, rT{order})
        end
    end

    # ordered ones should have defaults so leave them in
    push!(rot_types, rT)
end


#########################################################################
# Check fixed relationships
#########################################################################

# Do no rotation
R = eye(RotMatrix{Float64})
verbose_rot_test ? println("********************************\nIndentity checks\n********************************\n") : nothing
for rT in rot_types
    verbose_rot_test ? println(rT) : nothing
    rot_var = rT(R)
    null_var = eye(rT)
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

    # export to mutable
    mvu = Vector(rot_var)
    @contents_approx_eq(rot_var, mvu)

    # import from mutable
    rot_mvu = rT(mvu)
    @types_approx_eq(rot_var, rot_mvu)

    # export to immutable
    ivu = Vec(rot_var)
    @contents_approx_eq(rot_var, ivu)

    # import from immutable
    rot_ivu = rT(ivu)
    @types_approx_eq(rot_var, rot_ivu)

    # test typed stuff
    for eT in eltypes

        verbose_rot_test ? println("$(rT): $(eT)") : nothing

        # export to mutable
        mvt = Vector{eT}(rot_var)
        mvc = convert(Vector{eT}, mvu)
        @types_approx_eq(mvt, mvc)

        # import from mutable
        rot_mvt = rT(mvt)
        @contents_approx_eq_notype(rot_var, rot_mvt)

        # export to immutable
        ivt = Vec{numel(rT), eT}(rot_var)
        ivc = convert(Vec{numel(rT), eT}, ivt)
        @types_approx_eq(ivt, ivc)

        # import from immutable
        rot_ivt = rT(ivt)
        @contents_approx_eq_notype(rot_var, rot_ivt)


        # and test the element conversion on the rotation parameterization directly
        if (length(rT.parameters) == 1) || (TypeVar != typeof(rT.parameters[end-1]))
            rot_c = convert(rT{eT}, rot_var)
            @contents_approx_eq(rot_c, rot_mvt)
        end
    end
end

#########################################################################
# Check that Quaternions can import a point
#########################################################################

verbose_rot_test ? println("********************************\nQuaternion point import checks\n********************************\n") : nothing
@types_approx_eq(Quaternion([1.0,2.0,3.0]), Quaternion(0.0, 1.0, 2.0, 3.0))
@types_approx_eq(Quaternion(Vec(1.0,2.0,3.0)), Quaternion(0.0, 1.0, 2.0, 3.0))


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
            X = convert(rT_in, q)

            # round trip conversion
            Xd = convert(rT_in, convert(rT_out, X))

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
        Xo_t = rotate(R, X)
        @types_approx_eq_eps(Xo_t, Xo, 1e-10)

    end
end

