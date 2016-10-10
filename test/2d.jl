using Rotations, StaticArrays, Base.Test

@testset "2d Rotations" begin

    ###############################
    # Check fixed relationships
    ###############################

    @testset "Identity rotation checks" begin
        I = eye(SMatrix{2,2,Float64})
        I32 = eye(SMatrix{2,2,Float32})
        R = RotMatrix{2}
        @test eye(R)::R == I
        @test eye(R{Float32})::R{Float32} == I32
    end

    ################################
    # check on the inverse function
    ################################

    @testset "Testing inverse()" begin
        repeats = 100
        R = RotMatrix{2,Float64}
        I = eye(R)
        srand(0)
        for i = 1:repeats
            r = rand(R)
            @test inv(r) == r'
            @test inv(r) == r.'
            @test inv(r)*r ≈ I
            @test r*inv(r) ≈ I
        end
    end

    #########################################################################
    # Rotate some stuff
    #########################################################################

    # a random rotation of a random point
    @testset "Rotate Points" begin
        repeats = 100
        R = RotMatrix{2}
        srand(0)
        for i = 1:repeats
            r = rand(R)
            m = SMatrix(r)
            v = randn(SVector{2})

            @test r*v ≈ m*v
        end

        # Test Base.Vector also
        r = rand(R)
        m = SMatrix(r)
        v = randn(2)

        @test r*v ≈ m*v
    end

    # compose two random rotations
    @testset "Compose rotations" begin
        repeats = 100
        R = RotMatrix{2}
        srand(0)
        for i = 1:repeats
            r1 = rand(R)
            m1 = SMatrix(r1)

            r2 = rand(R)
            m2 = SMatrix(r2)

            @test r1*r2 ≈ m1*m2
        end
    end


    #########################################################################
    # Test conversions between rotation types
    #########################################################################
    @testset "Convert rotations" begin
        repeats = 100
        R = RotMatrix{2}
        srand(0)
        r1 = rand(R)
        @test R(r1) == r1
    end

    @testset "Types and products" begin
        for (R,T) in ((RotMatrix(pi/4), Float64),
                      (RotMatrix(Float32(pi/4)), Float32),
                      (RotMatrix{2,Float32}(pi/4), Float32))
            @test eltype(R) == T
            @test size(R) == (2,2)
            @test R == T[cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)]
            @test R * R ≈ T[0 -1; 1 0]
        end
    end

end

nothing
