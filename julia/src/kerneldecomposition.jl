#=
kerneldecomposition:
- Julia version: 
- Author: qd4314
- Date: 2019-03-19
=#

using KitBase

function metacheck(N)
    for n = 1:N
        checkexpansionHG(n)
    end
end

function checkexpansionHG(n, g = 0.3)
    ntheta = 1000
    theta = range(0, stop = 2 * pi, length = ntheta)
    hg(g, mu) = 1 / (4 * pi) * (1 - g^2) ./ (1 + g^2 - 2 * g * mu) .^ (3 / 2)
    y = hg.(g, cos.(theta))
    hgcoeffs = [(2 * l + 1) / (4 * pi) * g^l for l = 0:n]
    v = sf_legendre_Pl_array.(n, cos.(theta))
    yapprox = [sum(dot(v[i], hgcoeffs)) for i = 1:ntheta]
    println(norm(yapprox .- y) / norm(y))
end

function metadecomposition(N)
    err = zeros(N)
    for n = 1:N
        _, _, _, err[n] = kerneldecomposition(n)
    end
    return (err)
end

function kerneldecomposition(KS::KineticSolver, order = 20)
    g = 0.01
    HG = computeScatteringKernelHenyeyGreenstein(KS, g)


    index(l, m) = l * (order + 1) + m + 1

    Ωs = KS.Q.pointsxyz
    nq = KS.nquadpoints
    weights = KS.Q.weights

    Norder = (order + 1) * (order + 1)
    S = zeros(Norder, Norder)
    counter = 1
    for l = 0:order
        for m = -l:+l
            S[counter, counter] = g^l
            counter += 1
        end
    end


    Y = zeros(order + 1, 2 * (order + 1) + 1, nq)
    YY = zeros(Norder, nq)
    @polyvar x y z
    counter = 1
    for l = 0:order
        for m = -l:l
            sphericalh = ylm(l, m, x, y, z)
            for q = 1:nq
                Y[l+1, m+l+1, q] = sphericalh(x => Ωs[q, 1], y => Ωs[q, 2], z => Ωs[q, 3])
                YY[counter, q] = sphericalh(x => Ωs[q, 1], y => Ωs[q, 2], z => Ωs[q, 3])
            end
            counter += 1
        end
    end
    normY = zeros(Norder)
    for i = 1:Norder
        for k = 1:nq
            normY[i] = normY[i] + YY[i, k]^2
        end
    end
    normY = sqrt.(normY)


    Mat = zeros(nq, nq)
    O = zeros(nq, Norder)
    M = zeros(Norder, nq)
    for i = 1:nq
        for j = 1:nq
            val = 0
            counter = 1
            for l = 0:order
                for m = -l:l
                    val += g^l * Y[l+1, m+l+1, i] * Y[l+1, m+l+1, j] * weights[j]
                    O[i, counter] = YY[counter, i] # / normY[counter] 
                    M[counter, j] = YY[counter, j] * weights[j]# / normY[counter]
                    counter += 1
                end
            end
            Mat[i, j] = val
        end
    end


    v = O * S * M * ones(nq)

    counter = 1
    for q = 1:nq
        O[q, :] /= v[q]
    end

    #println(maximum(abs.(v.-HG*ones(nq))))

    mystring = string(norm((O * S * (M) .- HG)) / norm(HG))
    println("Error for HG - Prod(Matrixdecomposition)= " * mystring)
    err = norm(O * S * M .- HG) / norm(HG)
    return O, S, M, err

end

function kerneldecomposition(order = 20)
    KS = KineticSolver("config.txt")
    kerneldecomposition(KS, order)
end
