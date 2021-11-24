#=
computeKernels:
- Julia version:
- Author: qd4314
- Date: 2019-03-21
=#

function computeScatteringKernelConvolution(KS::KineticSolver)
    Nq = KS.nquadpoints
    weightMatrix = zeros(Nq, Nq)
    for q = 1:Nq
        Omega = KS.Q.pointsxyz[q, :]
        for l = 1:Nq
            Omegaprime = KS.Q.pointsxyz[l, :]
            #eps = KS.convolutionwidth/Nq^2 # do we need a second parameter to adjust how isotropic the scattering is?
            eps = KS.convolutionwidth / Nq
            #x = getdist(Omega,Omegaprime,"sphere")
            x = (1 - dot(Omega, Omegaprime))
            val = (KS.Q.weights[l] * 1.0 / eps * exp(-x^2 / (eps^2)))
            if (val < 1e-8)
                weightMatrix[q, l] = 0.0
            else
                weightMatrix[q, l] = val
            end
        end
        #weightMatrix[q, q] -= sum(weightMatrix[q, :])
        #weightMatrix[q,:] = weightMatrix[q,:]./KS.Q.weights[:]
        weightMatrix[q, :] = weightMatrix[q, :] ./ (sum(weightMatrix[q, :])) # normalize
        weightMatrix[q, :] = weightMatrix[q, :] ./ KS.Q.weights[:] # ensure that weightMatrix*KS.Q.weights[:] = ones(nq)
    end
    weightMatrix = sparse(weightMatrix)

    KS.OutscatteringVector = zeros(KS.nquadpoints)
    for k = 1:KS.nquadpoints
        KS.OutscatteringVector[k] = sum(weightMatrix[k, :])
    end

    return weightMatrix
end

function numericnnz(A, treshold = 1e-16)
    # Returns the number of nonzero entries of 'A'
    # A nonzero entry is any entry that has an
    # absolute value above the specified 'treshold'.
    count = 0
    for i = 1:prod(size(A))
        count += (abs(A[i]) < treshold)
    end
    return count
end


function analyze_nnz_for_convwidth(KS)
    N = zeros(0)
    I = zeros(0)
    nelements = 0
    for i = 10:10:1000
        KS.convolutionwidth = i
        S = computeScatteringKernelConvolution(KS)
        nelements = prod(size(S))
        n = nnz(S, 1e-8)
        append!(N, n)
        append!(I, i)
        println("$i, $n")
    end
    lineplot(I, N / nelements, xlabel = "convolutionwidth", ylabel = "ratio of nnz")
    return I, N
end


function computeScatteringKernelHenyeyGreenstein(KS::KineticSolver, g = 0.00000001)
    Nq = KS.nquadpoints
    HGScattering = zeros(Nq, Nq)
    hg(g, mu) = 1 / 4 / pi * (1 - g^2) ./ (1 + g^2 - 2 * g * mu) .^ (3 / 2)
    for q = 1:Nq
        for l = 1:Nq
            Omega = KS.Q.pointsxyz[q, :]
            Omegaprime = KS.Q.pointsxyz[l, :]
            mu = dot(Omega, Omegaprime)

            val = hg(g, mu)
            HGScattering[q, l] = ifelse(val < 1e-12, 0.0, val) * KS.Q.weights[l]
        end
        #	HGScattering[q, :] /= sum(HGScattering[q, :])##
        #	HGScattering[q,:] .* = KS.Q.weights[:] / sum(KS.Q.weights[:])
        #HGScattering[q,q] -= 1
    end
    #print((HGScattering*ones(Nq))
    v = HGScattering * ones(Nq)
    for q = 1:Nq
        for qprime = 1:Nq
            HGScattering[q, qprime] /= v[q]
        end
    end

    return HGScattering
end
