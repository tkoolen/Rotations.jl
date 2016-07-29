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


    #########################################################################
    # Test conversions between rotation types
    #########################################################################
    @testset "Convert rotations" begin
        repeats = 100
        @testset "convert $(R1) -> $(R2)" for R1 in all_types, R2 in rot_types
            srand(0)
            for i = 1:repeats
                r1 = rand(R1)
                m1 = SMatrix(r1)

                r2 = R2(r1)

                @test r2 ≈ m1
            end
        end
    end


    #########################################################################
    # Check angle and axis and inv work as expected
    #########################################################################

    @testset "Testing angle / axis extraction" begin
        repeats = 100
        @testset "$(R)" for R in rot_types
            srand(0)
            for i = 1:repeats
                r1 = rand(AngleAxis)

                angle = rotation_angle(r1)
                axis = rotation_axis(r1)

                r2 = R(r1)

                @test rotation_angle(r2) ≈ angle
                @test rotation_axis(r2) ≈ axis
            end
        end

        # TODO RotX, RotXY?
    end


    #########################################################################
    # Check roll, pitch, yaw constructors
    #########################################################################

    @testset "Testing roll / pitch / yaw constructors" begin
        repeats = 100
        s = 1e-4
        @testset "$(R)" for R in taitbyran_types
            srand(0)
            for i = 1:repeats
                roll = s*(rand()-0.5)
                pitch = s*(rand()-0.5)
                yaw = s*(rand()-0.5)

                # This tests whether the rotations are the same to first order
                # in roll, pitch and yaw. Second-order terms are small enough
                # to pass the isapprox() test, but not first-order terms.
                r1 = RotXYZ(roll, pitch, yaw)
                r2 = R(roll=roll, pitch=pitch, yaw=yaw)

                @test r1 ≈ r2
            end
        end
    end
end
