module Domaca01
using LinearAlgebra, Plots
import Base: \, adjoint, transpose

"""
    A = (V, I)

Za razpršene matrike definiramo novo strukturo RazprsenaMatrika, ki je hranjena v matrikah `V` in `I`.
Matrika zaradi prostorskih zahtev hrani vrednosti v matrikah `V` in `I`, velikosti n x m. Pri tem velja:
`V(i,j) = A(i, I(i, j))`.
"""
struct RazprsenaMatrika
    V
    I
end

"""
    ix = getIndex(A, i, j)

Vrne vrednost, ki se v matriki `A` nahaja na mestu `i, j`.
"""
function getIndex(A::RazprsenaMatrika, i, j)
    idx = findfirst(isequal(j), A.I[i])
    return idx === nothing ? 0.0 : A.V[i][idx]
end
Base.getindex(A::RazprsenaMatrika, i, j) = getIndex(A, i, j)

"""
    ix = setIndex(A, x, i, j)

Na mesto `i, j` v matriki `A` vnese vrednost `x`.
"""
function setIndex(A::RazprsenaMatrika, i, j, e)
    idx = findfirst(isequal(j), A.I[i])
    if idx === nothing
        ix = searchsortedlast(A.I[i], j)
        insert!(A.I[i], ix + 1, j)
        insert!(A.V[i], ix + 1, e)
    else
        A.V[i][idx] = e
    end
end
Base.setindex!(A::RazprsenaMatrika, e, i, j) = setIndex(A, i, j, e)

"""
    first = firstIndex(A, i)

Vrne prvo vrednost v matriki `A`.
"""
function firstIndex(A::RazprsenaMatrika, i)
    return A.I[i][1]
end
Base.firstindex(A::RazprsenaMatrika, i) = firstIndex(A, i)


"""
    first = lastIndex(A, i)

Vrne zanjo vrednost v matriki `A`.
"""
function lastIndex(A::RazprsenaMatrika, i)
    return A.I[i][end]
end
Base.lastindex(A::RazprsenaMatrika, i) = lastIndex(A, i)

"""
    b = A  * x

Izvede množenje matrike `A` in vektorja `x` iz desne.
"""
function Base.:*(A::RazprsenaMatrika, x)
    n = length(A.V)
    b = zeros(Float64, n)

    for i in 1:n
        for (j, idx) in enumerate(A.I[i])
            b[i] += A.V[i][j] * x[idx]
        end
    end
    return b
end


"""
    x, iter = SOR(A, b, x0, omega, tol)

Izračunaj vrednost `x` in število iteracij `iter` za matriko A s SOR iteracijo.
Pri tem podamo vektor `b`, začetni približek `x0`, parameter SOR iteracije `omega` in
pogoj za zaključek `tol`. Nekatere vrednosti omege ne konvergirajo vedno, zato je število
iteracija omejeno na 1000.
"""
function SOR(A::RazprsenaMatrika, b, x0, omega, tol=1e-10)
    n = length(b)
    x = copy(x0)

    iter = 0

    while true
        x_old = copy(x)
        for i in 1:n
            sum1 = 0.0
            sum2 = 0.0

            for j in A.I[i]
                if j < i
                    sum1 += A[i, j] * x[j]
                elseif j > i
                    sum2 += A[i, j] * x_old[j]
                end
            end

            x[i] = (1 - omega) * x_old[i] + (omega / A[i, i]) * (b[i] - sum1 - sum2)
        end

        iter += 1

        if maximum(abs.(A * x - b)) <= tol
            break
        end

        if iter >= 1000
            break
        end
    end

    return x, iter
end

"""
    A = ustvariRazsprsenoMatriko(vozlisca, robovi, vrednosti)

Ustvari in vrne razpršeno matriko `A` iz podanih vozlišč, povezav in vrednosti
grafa.
"""
function ustvariRazsprsenoMatriko(vozlisca, robovi, vrednosti)
    n = length(vozlisca)
    V = Vector{Vector{Float64}}(undef, n)
    I = Vector{Vector{Int}}(undef, n)
    
    for i in 1:n
        V[i] = Float64[]
        I[i] = Int[]
    end
    
    m = length(robovi)
    for k in 1:m
        i, j = robovi[k]
        v = vrednosti[k]
        push!(V[i], v)
        push!(I[i], j)
    end
    
    return RazprsenaMatrika(V, I)
end

"""
    grafOmeg, minStIteracij, optimalniOmega, minX = optimalnaOmega(A, b, x0, stGrafa)

Poišče in vrne optimalno vrednost `omega`, ki vodi do najhitrejše konvergence SOR 
iteracije. Prav tako vrne tudi graf, število iteracij in rezultat SOR iteracije 
pri optimalni omegi.
Pri tem izvede SOR iteracijo za 20 vrednosti med 0.0 in 2.0.
"""
function optimalnaOmega(A::RazprsenaMatrika, b, x0, stGrafa)
    omege = range(0, stop=2.0, length=20)
    stIteracij = Float64[]
    xValues = Vector{Vector{Float64}}()

    for omega in omege
        x, iter = SOR(A, b, x0, omega)
        push!(stIteracij, iter)
        push!(xValues, x)
    end

    minIteracij, minIx = findmin(stIteracij)
    minOmega = omege[minIx]
    minX = xValues[minIx]

    if (stGrafa == 1)
        return scatter(omege, stIteracij, xlabel="Omega", ylabel="Število iteracij",
                       label="Primer $stGrafa", title="Odvisnost hitrosti konvergence od vrednosti omege"), minIteracij, minOmega, minX
    end
    return scatter!(omege, stIteracij, label="Primer $stGrafa"), minIteracij, minOmega, minX
end

"""
    Lred, bx, by = sistemLaplace(vozlišča, robovi, xy;
                                    fiksni = Dict{Int,Tuple{Float64,Float64}}())

Zgradi reducirano Laplaceovo matriko in dva desna člena za vgrajevanje grafa.
* `vozlišča` - vektor oznak vozlišč  
* `robovi`   - vektor parov (i, j) (neusmerjeni, dovoljeni so dvojniki)  
* `xy`       - slovar Dict{Int,Tuple{Float64,Float64}} z začetnimi koordinatami
* `fiksni`   - fiksne koordinate
"""
function sistemLaplace(vozlišča, robovi, xy; fiksni = Dict())
    n  = length(vozlišča)
    idx = Dict(v=>k for (k,v) in pairs(vozlišča))

    # degree and adjacency
    deg = zeros(Int, n)
    adj = [Int[] for _ in 1:n]
    for (u,v) in robovi
        iu, iv = idx[u], idx[v]
        push!(adj[iu], iv); push!(adj[iv], iu)
        deg[iu] += 1;     deg[iv] += 1
    end

    V = Vector{Vector{Float64}}(undef, n)
    I = Vector{Vector{Int}}(undef, n)
    for i in 1:n
        V[i] = Float64[]
        I[i] = Int[]
        push!(V[i], -deg[i]);   push!(I[i], i)
        for j in adj[i]
            push!(V[i], 1.0);   push!(I[i], j)
        end
    end
    L = RazprsenaMatrika(V, I)

    free   = setdiff(vozlišča, keys(fiksni))
    m      = length(free)
    Vf, If = Vector{Vector{Float64}}(undef, m), Vector{Vector{Int}}(undef, m)
    bx     = zeros(Float64, m)
    by     = zeros(Float64, m)

    for (row, v) in enumerate(free)
        i = idx[v]
        Vf[row] = Float64[]
        If[row] = Int[]
        for (val, j) in zip(L.V[i], L.I[i])
            if haskey(fiksni, vozlišča[j])
                bx[row] -= val * fiksni[vozlišča[j]][1]
                by[row] -= val * fiksni[vozlišča[j]][2]
            else
                push!(Vf[row], val)
                push!(If[row], findfirst(==(vozlišča[j]), free))
            end
        end
    end
    return RazprsenaMatrika(Vf, If), bx, by
end

"""
    backslash(A::RazprsenaMatrika, b::Vector{Float64})

Reši sistem enačb Ax = b, kjer je A razpršena matrika (tipa RazprsenaMatrika),
z uporabo metode SOR z začetnim približkom ničelnega vektorja in privzeto
relaksacijsko vrednostjo omega = 1.0.

Ta funkcija omogoča uporabo operatorja backslash tudi za lastni podatkovni tip matrike.
"""
function \(A::RazprsenaMatrika, b::Vector{Float64})
    x0 = zeros(length(b))
    x, _ = SOR(A, b, x0, 1.0)
    return x
end

"""
    adjoint(A::RazprsenaMatrika)

Vrne kompleksno konjugirano transponirano matriko A.
Ker matrika vsebuje realne vrednosti, je to enakovredno transponiranju.
"""
function adjoint(A::RazprsenaMatrika)
    transpose(A)
end

"""
    transpose(A::RazprsenaMatrika)

Vrne transponirano matriko tipa RazprsenaMatrika. To pomeni, da se vrednosti iz vrstic
preslikajo v ustrezne stolpce nove matrike. Vrstice in stolpci se zamenjajo.
"""
function transpose(A::RazprsenaMatrika)
    nrows = length(A.V)
    all_cols = unique(vcat(A.I...))
    ncols = maximum(all_cols)

    Vt = [Float64[] for _ in 1:ncols]
    It = [Int[] for _ in 1:ncols]

    for i in 1:nrows
        for (val, j) in zip(A.V[i], A.I[i])
            push!(Vt[j], val)
            push!(It[j], i)
        end
    end

    return RazprsenaMatrika(Vt, It)
end

export RazprsenaMatrika, getIndex, setIndex, firstIndex, lastIndex, SOR, ustvariRazsprsenoMatriko, optimalnaOmega, sistemLaplace

end # module Domaca01
