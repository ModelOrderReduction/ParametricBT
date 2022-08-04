include("../src/ParametricBT.jl")
using .ParametricBT
using Random
using Test

function generate_test_system(n, m, p)
    A0 = rand(MersenneTwister(0), n, n);
    A0 = -A0'*A0
    A1 = rand(MersenneTwister(1), n, n);
    A1 = -A1'*A1
    A2 = rand(MersenneTwister(2), n, n);
    A2 = -A2'*A2
    A(x) = A0+A1*x+A2*x^2
    B(x) = rand(MersenneTwister(3), n, m)
    C(x) = rand(MersenneTwister(3), p, n)
    return A, B, C
end


@testset "ParametricBT.jl" begin
    f(x) = x^3
    A0 = rand(MersenneTwister(0), 5, 5);
    A1 = rand(MersenneTwister(1), 5, 5);
    A2 = rand(MersenneTwister(2), 5, 5);
    A(x) = A0+A1*x+A2*x^2
    x0 = 0.0
    Acoeffs = ParametricBT.compute_taylor_coeffs(A, x0, 2)
    @test length(Acoeffs) == 3
    @test Acoeffs[1] ≈ A0
    @test Acoeffs[2] ≈ A1
    @test Acoeffs[3] ≈ A2
end

    A, B, C = generate_test_system(5, 2, 3)
    ParametricBT.balancing(A, B, C, 2, 0.0)
