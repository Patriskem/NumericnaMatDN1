using Test, Domaca02, QuadGK, Distributions

tol = 1e-10

function cdfCalc(x)
    dist = Normal(0, 1)
    res = cdf(dist, x)
    return res
end

@testset "Gostota porazdelitve standardne normalne slučajne spremenljivke" begin
    # Primerjava izračunanih vrednosti z znanimi
    @test abs(normalPorGost(0.0) - 0.3989422804) <= tol
    @test abs(normalPorGost(1.0) - 0.2419707245) <= tol
    @test abs(normalPorGost(-1.0) - 0.2419707245) <= tol
    @test abs(normalPorGost(2.0) - 0.0539909665) <= tol
    @test abs(normalPorGost(-2.0) - 0.0539909665) <= tol

    # Primerjava integrala gostotne funkcije (mora biti 1)
    res, _ = quadgk(normalPorGost, -Inf, Inf)
    @test abs(res - 1) <= tol
end

@testset "Simpsonovo pravilo" begin
    # Integral funkcije f(x) = 1 na intervalu [a, b]
    f1(x) = 1
    a1, b1 = 0, 1
    pricakovan1 = b1 - a1
    res1 = simpsonovoPravilo(f1, a1, b1, 100)
    @test isapprox(res1, pricakovan1)

    # Integral funkcije f(x) = x na intervalu [0, 1]
    f2(x) = x
    a2, b2 = 0, 1
    pricakovan2 = 0.5
    res2 = simpsonovoPravilo(f2, a2, b2, 100)
    @test isapprox(res2, pricakovan2)

    # Integral funkcije f(x) = sin(x) na intervalu [0, π]
    f3(x) = sin(x)
    a3, b3 = 0, π
    pricakovan3 = 2  # Analitična vrednost integrala f(x) = sin(x) na [0, π] je 2
    res3 = simpsonovoPravilo(f3, a3, b3, 100)
    @test isapprox(res3, pricakovan3)
end

@testset "Normalna porazdelitev" begin
    @test abs(cdfCalc(-3.0) - normalPor(-3.0)) <= tol
    @test abs(cdfCalc(-2.0) - normalPor(-2.0)) <= tol
    @test abs(cdfCalc(-1.0) - normalPor(-1.0)) <= tol
    @test abs(cdfCalc(0.0) - normalPor(0.0))  <= tol
    @test abs(cdfCalc(1.0) - normalPor(1.0)) <= tol
    @test abs(cdfCalc(2.0) - normalPor(2.0)) <= tol
    @test abs(cdfCalc(3.0) - normalPor(3.0)) <= tol
end

@testset "Gauss-Legendre 2 točki" begin
    f(x) = x^2
    a, b = 0.0, 2.0
    pricakovan = (b^3 - a^3) / 3.0
    @test abs(gaussLegendre2P(f, a, b) - pricakovan) <= tol
end

@testset "Sestavljen Gauss-Legendre" begin
    f(x) = x^2
    a, b = 0.0, 2.0
    n = 4
    pricakovan = (b^3 - a^3) / 3.0
    @test abs(gaussLegendre2P(f, a, b) - pricakovan) <= tol
end

@testset "Primer integrala" begin
    f(x) = sin(x) / x
    a, b = 0.0, 5.0
    n = 10
    tol = 1e-10
    res, _ = primerGL()
    pricakovan, _ = quadgk(f, a, b, rtol=tol)
    @test isapprox(res, pricakovan)
end
