# function to perform tests of the rotation functions in the Rotations module
using StaticArrays
using Rotations

#=
########################
# Define helper methods
########################

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
=#

#####################################################################################
# build a full list of rotation types including the different angle ordering schemas
#####################################################################################

rot_types = (RotMatrix, Quat, SPQuat, AngleAxis, RodriguesVec,
             RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
             RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ)

one_types = (RotX, RotY, RotZ)
two_types = (RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)
taitbyran_types = (RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX)
all_types = (RotMatrix, Quat, SPQuat, AngleAxis, RodriguesVec,
             RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
             RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
             RotX, RotY, RotZ,
             RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)

###############################
# Start testing
###############################

@testset "Rotations Tests" begin

    ###############################
    # Check fixed relationships
    ###############################

    @testset "Identity rotation checks" begin
        I = eye(SMatrix{3,3,Float64})
        I32 = eye(SMatrix{3,3,Float32})
        @testset "$(R)" for R in all_types
            @test eye(R)::R == I
            @test eye(R{Float32})::R{Float32} == I32
        end
    end

    ################################
    # check on the inverse function
    ################################

    @testset "Testing inverse()" begin
        repeats = 100
        I = eye(RotMatrix{Float64})
        @testset "$(R)" for R in all_types
            srand(0)
            for i = 1:repeats
                r = rand(R)
                @test inv(r) == r'
                @test inv(r) == r.'
                @test inv(r)*r ≈ I
                @test r*inv(r) ≈ I
            end
        end
    end

    #########################################################################
    # Rotate some stuff
    #########################################################################

    # a random rotation of a random point
    @testset "Rotate Points" begin
        repeats = 100
        @testset "$(R)" for R in all_types
            srand(0)
            for i = 1:repeats
                r = rand(R)
                m = SMatrix(r)
                v = randn(SVector{3})

                @test r*v ≈ m*v
            end

            # Test Base.Vector also
            r = rand(R)
            m = SMatrix(r)
            v = randn(3)

            @test r*v ≈ m*v
        end
    end

    # compose two random rotations
    @testset "Compose rotations" begin
        repeats = 100
        @testset "$(R1) * $(R2)" for R1 in all_types, R2 in all_types
            srand(0)
            for i = 1:repeats
                r1 = rand(R1)
                m1 = SMatrix(r1)

                r2 = rand(R2)
                m2 = SMatrix(r2)

                @test r1*r2 ≈ m1*m2
            end
        end
    end


#=
    #########################################################################
    # Test conversions between rotation types
    #########################################################################

    # test random round trip conversions
    @testset "Rotation parameterization conversion checks" begin
        repeats, thresh = 100, 1e-12
        I = eye(RotMatrix{Float64})
        for rT_in in rot_types
            srand(0)

            # straight conversions
            @testset "$(rT_in) -> $(rT_out)" for rT_out in rot_types
                for i = 1:repeats                               # and each test

                    # start with a random quaternion
                    q = i == 1 ? Quaternion(1.0, 0.0, 0.0, 0.0) : nquatrand()  # the identity is a bit special for most parametrizations
                    X = convert(rT_in, q)

                    # round trip conversion
                    Xd = convert(rT_in, convert(rT_out, X))

                    # compare rotations before and after the round trip
                    Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the identity
                    rd = vecnorm(I - Rout)
                    @test rd <= thresh
                end
            end

            # element type conversion
            if (length(rT_in.parameters) == 1) || (TypeVar != typeof(rT_in.parameters[end-1]))          # dont try when there'a missing order parameter
                @testset "$(rT_in) -> $(rT_out) with element type conversion" for rT_out in rot_types[!ordered_type]

                    # start with a random quaternion
                    q = Quaternion(1.0, 0.0, 0.0, 0.0) # the identity is a bit special for most parametrizations
                    X = convert(rT_in{Float64}, q)

                    # check type inference without any element type conversions
                    @inferred convert(rT_out, X)
                    @inferred convert(rT_in, convert(rT_out, X))

                    # check type inference with element type conversions
                    @inferred convert(rT_out{Float32}, X)
                    @inferred convert(rT_in{Float64}, convert(rT_out{Float32}, X))
                    Xd = convert(rT_in{Float64}, convert(rT_out{Float32}, X))

                    # compare rotations before and after the round trip
                    Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the identity
                    rd = vecnorm(I - Rout)
                    @test rd <= thresh
                end
            end
        end
    end

    =#

    #=
    #########################################################################
    # Check angle and axis and inv work as expected
    #########################################################################

    @testset "Testing angle / axis extraction" begin
        theta_v, avec = linspace(-2*pi, 2*pi, 17), Vec(1.0, 0.0, 0.0)
        thresh = 1e-9
        @testset "$(rT)" for rT in [AngleAxis, RodriguesVec, Quaternion, SpQuat]
            for theta in theta_v
                aa = AngleAxis(theta, avec...)
                rot_var = rT(aa)
                theta_ex, avec_ex = rotation_angle(rot_var), rotation_axis(rot_var)

                # N.B. its OK of the axis is reversed and theta = 2*pi - theta
                dp = dot(avec_ex, avec)
                # @test_approx_eq_eps abs(dp)  1 thresh
                if (dp > 0)
                   @test abs(Rotations.wrap_angle(theta_ex - theta)) <= thresh
                else
                   @test abs(Rotations.wrap_angle((2*pi - theta_ex) - theta)) <= thresh
                end
            end
        end

        # make sure a slightly unnormalized Quaternion doesn't throw when extracting the angle
        q = Quaternion(1.0000000000000002, -5.040577330528428e-13, 5.663457441733299e-13, 2.452760217153127e-13)
        @test abs(rotation_angle(q)) < thresh

    end
    =#
end
