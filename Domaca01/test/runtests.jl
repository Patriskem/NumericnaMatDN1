using Test, Domaca01, LinearAlgebra

V = [[1, 2], [3, 4, 5]]
I = [[1, 2], [1, 2, 3]]
m = Domaca01.RazprsenaMatrika(V, I)

@testset "Razpršena matrika" begin
    @test isapprox(m.V, V)
    @test isapprox(m.I, I)
end

@testset "getIndex" begin
    @test Domaca01.getIndex(m, 1, 2) == 2
    @test Domaca01.getIndex(m, 1, 3) == 0
    @test Domaca01.getIndex(m, 2, 1) == 3
end

@testset "setIndex" begin
    Domaca01.setIndex(m, 1, 3, 6)
    Domaca01.setIndex(m, 1, 2, -2)
    Domaca01.setIndex(m, 2, 1, 0)
    @test Domaca01.getIndex(m, 1, 3) == 6
    @test Domaca01.getIndex(m, 1, 2) == -2
    @test Domaca01.getIndex(m, 2, 1) == 0
end

@testset "Prvi in zadnji indeks" begin
    @test Domaca01.firstIndex(m, 1) == 1
    @test Domaca01.lastIndex(m, 2) == 3
end

@testset "Množenje z desne z vektorjem" begin
    V = [[1, 2], [3, 4, 5]]
    I = [[1, 2], [1, 2, 3]]
    m = Domaca01.RazprsenaMatrika(V, I)
    x = [1.0, 2.0, 3.0]
    x_pricakovan = [1 * 1 + 2 * 2, 3 * 1 + 4 * 2 + 5 * 3]
    @test isapprox(m * x, x_pricakovan)
end

@testset "SOR iteracija" begin
    V = [[2.0], [1.0, 1.0], [3.0]]
    I = [[1], [2, 3], [3]]
    m = Domaca01.RazprsenaMatrika(V, I)
    b = [1.0, 4.0, 2.0]
    x0 = [0.0, 0.0, 0.0]
    omega = 1.25

    # Izračunamo pričakovan rezultat s pomočjo LU razcepa
    n = length(m.V)
    dense = zeros(Float64, n, n)
    for i in 1:n
        for (val, idx) in zip(m.V[i], m.I[i])
            dense[i, idx] = val
        end
    end
    F = lu(dense)

    x_pricakovan = F \ b
    x, iter = Domaca01.SOR(m, b, x0, omega)
    @test isapprox(x, x_pricakovan)
end

@testset "Pretvorba graf -> razpršena matrika" begin
    vozlisca = [1, 2, 3]
    robovi = [(1, 1), (2, 2), (2, 3), (3, 3)]
    vrednosti = [2.0, 1.0, 1.0, 3.0]
    A = ustvariRazsprsenoMatriko(vozlisca, robovi, vrednosti)

    V = [[2.0], [1.0, 1.0], [3.0]]
    I = [[1], [2, 3], [3]]
    @test isapprox(A.V, V)
    @test isapprox(A.I, I)
end

@testset "Laplace sistem" begin
    vozlišča = [1,2,3,4]
    robovi   = [(1,2),(2,3),(3,4)]
    fiksni   = Dict(1 => (-4.0,-5.0),
                    4 => ( 2.0, 3.0))
    L, bx, by = Domaca01.sistemLaplace(vozlišča, robovi, Dict(); fiksni = fiksni)
    @test length(L.V) == 2
    @test isapprox(bx, [ 4.0, -2.0])
    @test isapprox(by, [ 5.0, -3.0])
end

@testset "Vlaganje grafa" begin
    # Pricakovan rezultat
    V = [[2.0], [1.0, 1.0], [3.0]]
    I = [[1], [2, 3], [3]]
    m = Domaca01.RazprsenaMatrika(V, I)
    b = [1.0, 4.0, 2.0]

    n = length(m.V)
    dense = zeros(Float64, n, n)
    for i in 1:n
        for (val, idx) in zip(m.V[i], m.I[i])
            dense[i, idx] = val
        end
    end
    F = lu(dense)
    x_pricakovan = F \ b

    # Dejanski rezultat
    vozlisca = [1, 2, 3]
    robovi = [(1, 1), (2, 2), (2, 3), (3, 3)]
    vrednosti = [2.0, 1.0, 1.0, 3.0]
    A = ustvariRazsprsenoMatriko(vozlisca, robovi, vrednosti)  
    b = [1.0, 4.0, 2.0]
    x0 = [0.0, 0.0, 0.0]
    plot, minIteracij, minOmega, rez = optimalnaOmega(A, b, x0, 25)  # povečaj število iteracij
    @test isapprox(rez, x_pricakovan; atol=1e-3)  # dodaš toleranco, če je potrebna
end