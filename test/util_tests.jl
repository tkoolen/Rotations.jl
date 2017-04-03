@testset "Util" begin
    @testset "Perpendicular vector" begin
        for i = 1 : 100
            vec = rand(SVector{3, Float64})
            perp = Rotations.perpendicular_vector(vec)
            @test norm(perp) >= maximum(abs.(vec))
            @test isapprox(dot(vec, perp), 0.; atol = 1e-10)
        end
    end
end
