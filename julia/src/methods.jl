function methods!(
    KS::KineticSolver,
    psi::Array{Float64,3},
    psiOld::Array{Float64,3},
    k::Int64,
)


    if (KS.ThisQuadratureType == "octa" || KS.ThisQuadratureType == "ico") &&
       KS.rotationmagnitude > 0 # rSN rotation step
        u = randn(3)
        u = u / norm(u)
        theta = KS.rotationmagnitude * KS.dt / KS.norder
        rotateandinterpolate!(u, theta, KS.Q, psiOld, psi, "forward") # writes into psiOld
        @inbounds for i = 1:length(psi)
            psi[i] = psiOld[i]
        end

    elseif KS.ThisQuadratureType == "tensorized" && KS.rotationmagnitude > 0 # Rotation just around the z-axis
        if k % 2 == 0
            theta = KS.rotationmagnitude * KS.dt / KS.norder       # and ensure it lies in Delta psi
            alpha = theta / (pi / KS.norder) % 1.0
            rotateOnly!([0.0; 0.0; 1.0], theta, KS.Q, "forward")
            psiInterp = zeros(Int(KS.norder / 2), 2 * KS.norder)
            for i = 1:size(psi, 2)
                for j = 1:size(psi, 3)
                    psiMat = reshape(psi[:, i, j], Int(KS.norder / 2), 2 * KS.norder)
                    for q = 1:(2*KS.norder-1)
                        psiInterp[:, q] =
                            (1 - alpha) * psiMat[:, q] + alpha * psiMat[:, q+1]
                    end
                    psiInterp[:, 2*KS.norder] =
                        (1 - alpha) * psiMat[:, 2*KS.norder] + alpha * psiMat[:, 1]
                    psi[:, i, j] = psiInterp[:]
                end
            end
        else
            theta = KS.rotationmagnitude * KS.dt / KS.norder       # and ensure it lies in Delta psi
            alpha = theta % (pi / KS.norder)
            rotateOnly!([0.0; 0.0; -1.0], theta, KS.Q, "forward")
            psiInterp = zeros(Int(KS.norder / 2), 2 * KS.norder)
            for i = 1:size(psi, 2)
                for j = 1:size(psi, 3)
                    psiMat = reshape(psi[:, i, j], Int(KS.norder / 2), 2 * KS.norder)
                    for q = 2:2*KS.norder
                        psiInterp[:, q] =
                            (1 - alpha) * psiMat[:, q] + alpha * psiMat[:, q-1]
                    end
                    psiInterp[:, 1] =
                        (1 - alpha) * psiMat[:, 1] + alpha * psiMat[:, 2*KS.norder]
                    psi[:, i, j] = psiInterp[:]
                end
            end
        end
    end

    if KS.convolutionmagnitude > 0 # Idea of Cory Hauck, perform convolution in the angular domain

        @inbounds @views for j = 1:size(psi, 3)
            for i = 1:size(psi, 2)
                if norm(psi[:, i, j]) > 1e-12
                    mul!(
                        KS.tmp,
                        KS.ScatteringKernelConvolution,
                        KS.Q.weights .* psi[:, i, j],
                    )
                    psi[:, i, j] .+=
                        KS.dt .* KS.convolutionmagnitude .*
                        (KS.tmp - KS.OutscatteringVector .* KS.Q.weights .* psi[:, i, j])
                end
            end
        end
    end
end



function solveFluxUpwind!(KS::KineticSolver, psi::Array{Float64,3}, flux::Array{Float64,3})
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
        s3 = psi[q, i, j]
        northflux = 7 / 4 * s3 - s2 + 1 / 4 * s1
        southflux = 7 / 4 * s2 - s1 + 1 / 4 * s0

        s0 = psi[q, i-3, j]
        s1 = psi[q, i-2, j]
        s2 = psi[q, i-1, j]
        s3 = psi[q, i, j]
        eastflux = 7 / 4 * s3 - s2 + 1 / 4 * s1
        westflux = 7 / 4 * s2 - s1 + 1 / 4 * s0

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end
    #PosNeg
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxPosNeg
        s0 = psi[q, i, j+3]
        s1 = psi[q, i, j+2]
        s2 = psi[q, i, j+1]
        s3 = psi[q, i, j]
        northflux = 7 / 4 * s2 - s1 + 1 / 4 * s0
        southflux = 7 / 4 * s3 - s2 + 1 / 4 * s1

        s0 = psi[q, i-3, j]
        s1 = psi[q, i-2, j]
        s2 = psi[q, i-1, j]
        s3 = psi[q, i, j]
        eastflux = 7 / 4 * s3 - s2 + 1 / 4 * s1
        westflux = 7 / 4 * s2 - s1 + 1 / 4 * s0

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end

    # NegPos
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxNegPos
        s0 = psi[q, i, j-3]
        s1 = psi[q, i, j-2]
        s2 = psi[q, i, j-1]
        s3 = psi[q, i, j]
        northflux = 7 / 4 * s3 - s2 + 1 / 4 * s1
        southflux = 7 / 4 * s2 - s1 + 1 / 4 * s0

        s0 = psi[q, i+3, j]
        s1 = psi[q, i+2, j]
        s2 = psi[q, i+1, j]
        s3 = psi[q, i, j]
        westflux = 7 / 4 * s3 - s2 + 1 / 4 * s1
        eastflux = 7 / 4 * s2 - s1 + 1 / 4 * s0

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end

    # NegNeg
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxNegNeg
        s0 = psi[q, i, j+3]
        s1 = psi[q, i, j+2]
        s2 = psi[q, i, j+1]
        s3 = psi[q, i, j]
        northflux = 7 / 4 * s2 - s1 + 1 / 4 * s0
        southflux = 7 / 4 * s3 - s2 + 1 / 4 * s1

        s0 = psi[q, i+3, j]
        s1 = psi[q, i+2, j]
        s2 = psi[q, i+1, j]
        s3 = psi[q, i, j]
        westflux = 7 / 4 * s3 - s2 + 1 / 4 * s1
        eastflux = 7 / 4 * s2 - s1 + 1 / 4 * s0

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end
end


function Euler!(
    KS::KineticSolver,
    psi::Array{Float64,3},
    flux::Array{Float64,3},
    dt::Float64,
)
    # computes the Euler step
    @inbounds for j in KS.rangex, i in KS.rangey
        integral = 0.0
        @inbounds for q in KS.rangequad
            integral += psi[q, i, j] .* KS.Q.weights[q]
        end
        integral = integral * KS.SigmaS[i, j] / 4.0 / pi
        sigmaT = KS.SigmaT[i, j]
        @inbounds for q in KS.rangequad
            psi[q, i, j] += dt .* (-psi[q, i, j] .* sigmaT - flux[q, i, j] + integral)
        end
    end
end

# the first minmod code is the fast version of the second minmod  below that is commented
@inline minmod(x::Float64, y::Float64) = ifelse(x < 0, clamp(y, x, 0.0), clamp(y, 0.0, x))
#@inline function minmod(x::Float64, y::Float64)
#    return sign(x) * max(0.0, min(abs(x),y*sign(x) ) );
#end

@inline function slopefit(left::Float64, center::Float64, right::Float64)
    tmp = minmod(0.5 * (right - left), 2.0 * (center - left))
    return minmod(2.0 * (right - center), tmp)
end

function computeMoments(KS::KineticSolver, psi::Array{Float64,3})
    rho = zeros(KS.ny, KS.nx)
    for q in KS.rangequad, j = 1:KS.nx, i = 1:KS.ny#KS.rangex, i=KS.rangey
        rho[i, j] += psi[q, i, j] * KS.Q.weights[q]
    end
    return rho
end



function forcePeriodicity!(KS::KineticSolver, psi::Array{Float64,3})
    if KS.periodicXflag
        for q = 1:KS.nquadpoints
            for i = 1:KS.ny
                psi[q, i, 1] = psi[q, i, end-3]
                psi[q, i, 2] = psi[q, i, end-2]
                psi[q, i, end-1] = psi[q, i, 3]
                psi[q, i, end] = psi[q, i, 4]
            end
        end
    end
    if KS.periodicYflag
        for q = 1:KS.nquadpoints
            for j = 1:KS.nx
                psi[q, 1, j] = psi[q, end-3, j]
                psi[q, 2, j] = psi[q, end-2, j]
                psi[q, end-1, j] = psi[q, 3, j]
                psi[q, end, j] = psi[q, 4, j]
            end
        end
    end
end

function copyFromToNoGhosts!(
    psiFrom::Array{Float64,3},
    psiTo::Array{Float64,2},
    KS::KineticSolver,
)
    for q in KS.rangequad
        k = 0
        for i = 3:KS.ny-2
            for j = 3:KS.nx-3
                k = k + 1
                psiTo[q, k] = psiFrom[q, i, j]
            end
        end
    end
end

function copyFromToWithGhosts!(
    psiFrom::Array{Float64,2},
    psiTo::Array{Float64,3},
    KS::KineticSolver,
)
    for q in KS.rangequad
        k = 0
        for i = 1:KS.ny-4
            for j = 1:KS.nx-4
                k = k + 1
                psiTo[q, i+2, j+2] = psiFrom[q, k]
            end
        end
    end
end

function solveFlux!(KS::KineticSolver, phi::Array{Float64,3}, flux::Array{Float64,3})
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
        s1 = phi[q, i, j-2]
        s2 = phi[q, i, j-1]
        s3 = phi[q, i, j]
        s4 = phi[q, i, j+1]
        northflux = s3 + 0.5 .* slopefit(s2, s3, s4)
        southflux = s2 + 0.5 .* slopefit(s1, s2, s3)

        s1 = phi[q, i-2, j]
        s2 = phi[q, i-1, j]
        s3 = phi[q, i, j]
        s4 = phi[q, i+1, j]
        eastflux = s3 + 0.5 .* slopefit(s2, s3, s4)
        westflux = s2 + 0.5 .* slopefit(s1, s2, s3)

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end
    #PosNeg
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxPosNeg
        s1 = phi[q, i, j-1]
        s2 = phi[q, i, j]
        s3 = phi[q, i, j+1]
        s4 = phi[q, i, j+2]
        northflux = s3 - 0.5 .* slopefit(s2, s3, s4)
        southflux = s2 - 0.5 .* slopefit(s1, s2, s3)

        s1 = phi[q, i-2, j]
        s2 = phi[q, i-1, j]
        s3 = phi[q, i, j]
        s4 = phi[q, i+1, j]
        eastflux = s3 + 0.5 .* slopefit(s2, s3, s4)
        westflux = s2 + 0.5 .* slopefit(s1, s2, s3)

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end

    # NegPos
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxNegPos
        s1 = phi[q, i, j-2]
        s2 = phi[q, i, j-1]
        s3 = phi[q, i, j]
        s4 = phi[q, i, j+1]
        northflux = s3 + 0.5 .* slopefit(s2, s3, s4)
        southflux = s2 + 0.5 .* slopefit(s1, s2, s3)

        s1 = phi[q, i-1, j]
        s2 = phi[q, i, j]
        s3 = phi[q, i+1, j]
        s4 = phi[q, i+2, j]
        eastflux = s3 - 0.5 .* slopefit(s2, s3, s4)
        westflux = s2 - 0.5 .* slopefit(s1, s2, s3)

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end

    # NegNeg
    @inbounds for j in KS.rangex, i in KS.rangey, q in idxNegNeg
        s1 = phi[q, i, j-1]
        s2 = phi[q, i, j]
        s3 = phi[q, i, j+1]
        s4 = phi[q, i, j+2]
        northflux = s3 - 0.5 .* slopefit(s2, s3, s4)
        southflux = s2 - 0.5 .* slopefit(s1, s2, s3)

        s1 = phi[q, i-1, j]
        s2 = phi[q, i, j]
        s3 = phi[q, i+1, j]
        s4 = phi[q, i+2, j]
        eastflux = s3 - 0.5 .* slopefit(s2, s3, s4)
        westflux = s2 - 0.5 .* slopefit(s1, s2, s3)

        flux[q, i, j] =
            KS.Q.pointsxyz[q, 1] ./ KS.dx .* (eastflux - westflux) +
            KS.Q.pointsxyz[q, 2] ./ KS.dy .* (northflux - southflux)
    end
end
