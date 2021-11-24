#=
sweeping:
- Julia version: 
- Author: qd4314
- Date: 2019-03-19
=#



function sweeping(
    KS::KineticSolver,
    psi::Array{Float64,3},
    Source::Array{Float64,3},
    ScatteringKernel::Array{Float64,2},
)
    rhs = zeros(KS.nquadpoints, KS.nx, KS.ny)
    io = open("mytestfile.txt", "a")
    #println("")
    #println("sweeping is called")
    # loop over pseudo-time
    for k = 1:KS.nSweepingSteps
        # setup rhs
        print(" ")
        println(k)
        @inbounds for j in KS.rangex, i in KS.rangey
            rhs[:, i, j] =
                Source[:, i, j] +
                KS.convolutionmagnitude .* KS.kernelTimesWeightsSparse * (psi[:, i, j])
        end

        # check residual
        psiOld = deepcopy(psi)
        # perform sweep
        solveFluxSweeping!(KS, psi, rhs)

        # check residual
        #println("-> Sweeping norm: ",norm(psi-psiOld)*KS.contractionEstimate)
        if norm(psi - psiOld) * KS.contractionEstimate < KS.epsSweeping
            writedlm(io, k)
            close(io)
            return psi
        end
    end
    return psi
end

function sweepingWithKernels(
    KS::KineticSolver,
    psi::Array{Float64,3},
    Source::Array{Float64,3},
    asKernel::Array{Float64,2},
    physKernel::Array{Float64,2},
)
    rhs = zeros(KS.nquadpoints, KS.nx, KS.ny)
    io = open("mytestfile.txt", "a")
    #println("")
    #println("sweeping is called")
    # loop over pseudo-time
    for k = 1:KS.nSweepingSteps
        # setup rhs
        print(" ")
        println(k)
        @inbounds for j in KS.rangex, i in KS.rangey
            rhs[:, i, j] =
                Source[:, i, j] +
                KS.SigmaS[i, j] .* physKernel * psi[:, i, j] +
                KS.convolutionmagnitude .* asKernel * psi[:, i, j]
        end

        # check residual
        psiOld = deepcopy(psi)
        # perform sweep
        solveFluxSweeping!(KS, psi, rhs)
        #solveFluxSweepingFirstOrder!(KS,psi,rhs);

        # check residual
        println(
            "-> Sweeping norm: ",
            norm(psi - psiOld),
            " ; constraction: ",
            KS.contractionEstimate,
        )
        if norm(psi - psiOld) * KS.contractionEstimate < KS.epsSweeping
            writedlm(io, k)
            close(io)
            return psi
        end
    end
    return psi
end

function solveFluxSweepingBiggerStencil!(
    KS::KineticSolver,
    psi::Array{Float64,3},
    rhs::Array{Float64,3},
)
    # computes the numerical flux over cell boundaries for each ordinate
    # for faster computation, we split the iteration over quadrature points
    # into four different blocks: North West, Nort East, Sout West, South East
    # this corresponds to the direction the ordinates point to
    idxPosPos = findall((KS.Q.pointsxyz[:, 1] .>= 0.0) .& (KS.Q.pointsxyz[:, 2] .>= 0.0))
    idxPosNeg = findall((KS.Q.pointsxyz[:, 1] .>= 0.0) .& (KS.Q.pointsxyz[:, 2] .< 0.0))
    idxNegPos = findall((KS.Q.pointsxyz[:, 1] .< 0.0) .& (KS.Q.pointsxyz[:, 2] .>= 0.0))
    idxNegNeg = findall((KS.Q.pointsxyz[:, 1] .< 0.0) .& (KS.Q.pointsxyz[:, 2] .< 0.0))


    # PosPos
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxPosPos
        s0 = psi[q, i, j-3]
        s1 = psi[q, i, j-2]
        s2 = psi[q, i, j-1]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-s2 + 1 / 4 * s1)
        fluxExplY =
            fluxExplY - KS.Q.pointsxyz[q, 2] ./ KS.dy .* ((7 / 4) * s2 - s1 + 1 / 4 * s0)
        fluxImplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (7 / 4)

        s0 = psi[q, i-3, j]
        s1 = psi[q, i-2, j]
        s2 = psi[q, i-1, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-s2 + 1 / 4 * s1)
        fluxExplX =
            fluxExplX - KS.Q.pointsxyz[q, 1] ./ KS.dx .* ((7 / 4) * s2 - s1 + 1 / 4 * s0)
        fluxImplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (7 / 4)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end
    #PosNeg
    @inbounds for j in Iterators.reverse(KS.rangex), i in KS.rangey, q in idxPosNeg
        s0 = psi[q, i, j+1]
        s1 = psi[q, i, j+2]
        s2 = psi[q, i, j+3]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* ((7 / 4) * s0 - s1 + (1 / 4) * s2)
        fluxExplY = fluxExplY - KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-s0 + (1 / 4) * s1)
        fluxImplY = -KS.Q.pointsxyz[q, 2] ./ KS.dy .* (7 / 4)

        s0 = psi[q, i-3, j]
        s1 = psi[q, i-2, j]
        s2 = psi[q, i-1, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-s2 + 1 / 4 * s1)
        fluxExplX =
            fluxExplX - KS.Q.pointsxyz[q, 1] ./ KS.dx .* ((7 / 4) * s2 - s1 + 1 / 4 * s0)
        fluxImplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (7 / 4)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end

    # NegPos
    @inbounds for j in KS.rangex, i in Iterators.reverse(KS.rangey), q in idxNegPos
        s0 = psi[q, i, j-3]
        s1 = psi[q, i, j-2]
        s2 = psi[q, i, j-1]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-s2 + 1 / 4 * s1)
        fluxExplY =
            fluxExplY - KS.Q.pointsxyz[q, 2] ./ KS.dy .* ((7 / 4) * s2 - s1 + 1 / 4 * s0)
        fluxImplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (7 / 4)

        s0 = psi[q, i+1, j]
        s1 = psi[q, i+2, j]
        s2 = psi[q, i+3, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* ((7 / 4) * s0 - s1 + (1 / 4) * s2)
        fluxExplX = fluxExplX - KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-s0 + 1 / 4 * s1)
        fluxImplX = -KS.Q.pointsxyz[q, 1] ./ KS.dx .* (7 / 4)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end

    # NegNeg
    @inbounds for j in Iterators.reverse(KS.rangex),
        i in Iterators.reverse(KS.rangey),
        q in idxNegNeg

        s0 = psi[q, i, j+1]
        s1 = psi[q, i, j+2]
        s2 = psi[q, i, j+3]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* ((7 / 4) * s0 - s1 + (1 / 4) * s2)
        fluxExplY = fluxExplY - KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-s0 + (1 / 4) * s1)
        fluxImplY = -KS.Q.pointsxyz[q, 2] ./ KS.dy .* (7 / 4)

        s0 = psi[q, i+1, j]
        s1 = psi[q, i+2, j]
        s2 = psi[q, i+3, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* ((7 / 4) * s0 - s1 + (1 / 4) * s2)
        fluxExplX = fluxExplX - KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-s0 + 1 / 4 * s1)
        fluxImplX = -KS.Q.pointsxyz[q, 1] ./ KS.dx .* (7 / 4)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end
end

function solveFluxSweeping!(KS::KineticSolver, psi::Array{Float64,3}, rhs::Array{Float64,3})
    # computes the numerical flux over cell boundaries for each ordinate
    # for faster computation, we split the iteration over quadrature points
    # into four different blocks: North West, Nort East, Sout West, South East
    # this corresponds to the direction the ordinates point to
    idxPosPos = findall((KS.Q.pointsxyz[:, 1] .>= 0.0) .& (KS.Q.pointsxyz[:, 2] .>= 0.0))
    idxPosNeg = findall((KS.Q.pointsxyz[:, 1] .>= 0.0) .& (KS.Q.pointsxyz[:, 2] .< 0.0))
    idxNegPos = findall((KS.Q.pointsxyz[:, 1] .< 0.0) .& (KS.Q.pointsxyz[:, 2] .>= 0.0))
    idxNegNeg = findall((KS.Q.pointsxyz[:, 1] .< 0.0) .& (KS.Q.pointsxyz[:, 2] .< 0.0))


    # PosPos
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxPosPos
        s2 = psi[q, i, j-2]
        s1 = psi[q, i, j-1]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-2 * s1 + 1 / 2 * s2)
        fluxImplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (3 / 2)

        s2 = psi[q, i-2, j]
        s1 = psi[q, i-1, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-2 * s1 + 1 / 2 * s2)
        fluxImplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (3 / 2)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end
    #PosNeg
    @inbounds for j in Iterators.reverse(KS.rangex), i in KS.rangey, q in idxPosNeg
        s1 = psi[q, i, j+1]
        s2 = psi[q, i, j+2]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (2 * s1 - (1 / 2) * s2)
        fluxImplY = -KS.Q.pointsxyz[q, 2] ./ KS.dy .* (3 / 2)

        s2 = psi[q, i-2, j]
        s1 = psi[q, i-1, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-2 * s1 + 1 / 2 * s2)
        fluxImplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (3 / 2)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end

    # NegPos
    @inbounds for j in KS.rangex, i in Iterators.reverse(KS.rangey), q in idxNegPos
        s2 = psi[q, i, j-2]
        s1 = psi[q, i, j-1]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-2 * s1 + 1 / 2 * s2)
        fluxImplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (3 / 2)

        s1 = psi[q, i+1, j]
        s2 = psi[q, i+2, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (2 * s1 - (1 / 2) * s2)
        fluxImplX = -KS.Q.pointsxyz[q, 1] ./ KS.dx .* (3 / 2)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end

    # NegNeg
    @inbounds for j in Iterators.reverse(KS.rangex),
        i in Iterators.reverse(KS.rangey),
        q in idxNegNeg

        s1 = psi[q, i, j+1]
        s2 = psi[q, i, j+2]

        fluxExplY = KS.Q.pointsxyz[q, 2] ./ KS.dy .* (2 * s1 - (1 / 2) * s2)
        fluxImplY = -KS.Q.pointsxyz[q, 2] ./ KS.dy .* (3 / 2)

        s1 = psi[q, i+1, j]
        s2 = psi[q, i+2, j]

        fluxExplX = KS.Q.pointsxyz[q, 1] ./ KS.dx .* (2 * s1 - (1 / 2) * s2)
        fluxImplX = -KS.Q.pointsxyz[q, 1] ./ KS.dx .* (3 / 2)

        fluxImplicit = fluxImplY + fluxImplX

        flux = fluxExplX + fluxExplY

        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImplicit))
    end
end

function solveFluxSweepingFirstOrder!(
    KS::KineticSolver,
    psi::Array{Float64,3},
    rhs::Array{Float64,3},
)
    # computes the numerical flux over cell boundaries for each ordinate
    # for faster computation, we split the iteration over quadrature points
    # into four different blocks: North West, Nort East, Sout West, South East
    # this corresponds to the direction the ordinates point to
    idxPosPos = findall((KS.Q.pointsxyz[:, 1] .>= 0.0) .& (KS.Q.pointsxyz[:, 2] .>= 0.0))
    idxPosNeg = findall((KS.Q.pointsxyz[:, 1] .>= 0.0) .& (KS.Q.pointsxyz[:, 2] .< 0.0))
    idxNegPos = findall((KS.Q.pointsxyz[:, 1] .< 0.0) .& (KS.Q.pointsxyz[:, 2] .>= 0.0))
    idxNegNeg = findall((KS.Q.pointsxyz[:, 1] .< 0.0) .& (KS.Q.pointsxyz[:, 2] .< 0.0))


    # PosPos
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxPosPos
        s1 = psi[q, i, j-2]
        s2 = psi[q, i, j-1]
        s3 = psi[q, i, j]
        s4 = psi[q, i, j+1]
        northflux = s3#+0.5 .*slopefit(s2,s3,s4)
        southflux = s2#+0.5 .*slopefit(s1,s2,s3)

        s1 = psi[q, i-2, j]
        s2 = psi[q, i-1, j]
        s3 = psi[q, i, j]
        s4 = psi[q, i+1, j]
        eastflux = s3#+0.5 .*slopefit(s2,s3,s4)
        westflux = s2#+0.5 .*slopefit(s1,s2,s3)

        fluxImpllicit = KS.Q.pointsxyz[q, 1] ./ KS.dx + KS.Q.pointsxyz[q, 2] ./ KS.dy

        flux =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-southflux)
        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImpllicit))
    end
    #PosNeg
    @inbounds for j in Iterators.reverse(KS.rangex), i in KS.rangey, q in idxPosNeg
        s1 = psi[q, i, j-1]
        s2 = psi[q, i, j]
        s3 = psi[q, i, j+1]
        s4 = psi[q, i, j+2]
        northflux = s3#-0.5 .* slopefit(s2,s3,s4)
        southflux = s2#-0.5 .*slopefit(s1,s2,s3)

        s1 = psi[q, i-2, j]
        s2 = psi[q, i-1, j]
        s3 = psi[q, i, j]
        s4 = psi[q, i+1, j]
        eastflux = s3#+0.5 .*slopefit(s2,s3,s4)
        westflux = s2#+0.5 .*slopefit(s1,s2,s3)

        fluxImpllicit =
            KS.Q.pointsxyz[q, 1] ./ KS.dx + KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-1)
        flux =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux)
        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImpllicit))
    end

    # NegPos
    @inbounds for j in KS.rangex, i in Iterators.reverse(KS.rangey), q in idxNegPos
        s1 = psi[q, i, j-2]
        s2 = psi[q, i, j-1]
        s3 = psi[q, i, j]
        s4 = psi[q, i, j+1]
        northflux = s3#+0.5 .*slopefit(s2,s3,s4)
        southflux = s2#+0.5 .*slopefit(s1,s2,s3)

        s1 = psi[q, i-1, j]
        s2 = psi[q, i, j]
        s3 = psi[q, i+1, j]
        s4 = psi[q, i+2, j]
        eastflux = s3#-0.5 .*slopefit(s2,s3,s4)
        westflux = s2#-0.5 .*slopefit(s1,s2,s3)

        fluxImpllicit =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-1) + KS.Q.pointsxyz[q, 2] ./ KS.dy
        flux =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-southflux)
        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImpllicit))
    end

    # NegNeg
    @inbounds for j in Iterators.reverse(KS.rangex),
        i in Iterators.reverse(KS.rangey),
        q in idxNegNeg

        s1 = psi[q, i, j-1]
        s2 = psi[q, i, j]
        s3 = psi[q, i, j+1]
        s4 = psi[q, i, j+2]
        northflux = s3#-0.5 .* slopefit(s2,s3,s4)
        southflux = s2#-0.5 .* slopefit(s1,s2,s3)

        s1 = psi[q, i-1, j]
        s2 = psi[q, i, j]
        s3 = psi[q, i+1, j]
        s4 = psi[q, i+2, j]
        eastflux = s3#-0.5 .* slopefit(s2,s3,s4)
        westflux = s2#-0.5 .* slopefit(s1,s2,s3)

        fluxImpllicit =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (-1) + KS.Q.pointsxyz[q, 2] ./ KS.dy .* (-1)
        flux =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux)
        psi[q, i, j] =
            KS.dt .* (-flux + rhs[q, i, j]) /
            (1 + KS.dt .* (KS.SigmaT[i, j] + fluxImpllicit))
    end
end
