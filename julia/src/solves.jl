#=
solves:
- Julia version: 
- Author: qd4314
- Date: 2019-03-21
=#

function linFunction(phiVector::AbstractVector, KS::KineticSolver)

    # reshape vector input to tensor
    phi = reshape(phiVector, (KS.DecompNumberElements, KS.nx, KS.ny))

    OSigma = KS.OSigma
    M = KS.DecompM

    # compute sweeping rhs which is O*Sigma*phi
    OSigmaPhi = zeros(KS.nquadpoints, KS.nx, KS.ny)
    @einsum OSigmaPhi[q, i, j] = OSigma[q, moments] * phi[moments, i, j]

    # multiply with sigmaS
    @inbounds for j in KS.rangex, i in KS.rangey
        OSigmaPhi[:, i, j] = OSigmaPhi[:, i, j] .* KS.SigmaS[i, j]
    end

    # perform sweeping
    Y = sweeping(KS, deepcopy(KS.psiC), OSigmaPhi, KS.ScatteringKernelConvolution)

    #Y= sweeping(KS,0.0*deepcopy(KS.psiC),OSigmaPhi,KS.ScatteringKernelConvolution);  

    KS.psiC .= Y

    # multiply with M
    MY = zeros(KS.DecompNumberElements, KS.nx, KS.ny)
    @einsum MY[m, i, j] = M[m, q] * Y[q, i, j]

    # substract from phi to get A*phi output
    Aphi = phi .- MY

    # reshape for vector output
    AphiVector = reshape(Aphi, prod((KS.DecompNumberElements, KS.nx, KS.ny)))

    return AphiVector

end

function corysolve!(KS::KineticSolver)

    ### DEBUGGING
    # turn off scattering and absorption
    #KS.SigmaS = zeros(size(KS.SigmaS));
    #KS.SigmaT = zeros(size(KS.SigmaT));
    #KS.convolutionmagnitude = 0.0;
    close("all")
    io = open("mytestfile.txt", "a")

    writedlm(io, "start!!!!!!!!")
    close(io)
    # set convolution magnitude

    # set max number of sweeping steps
    KS.nSweepingSteps = 1000 # 1 is one-shot
    KS.epsSweeping = 1e-3#1e-7;

    # update SigmaT if convolution is active

    if KS.ictype == "linesource"
        println("Problem Type is ", KS.ictype)
        sigmaA = 0.0 # hard coded for linesource
    else
        sigmaA = KS.SigmaT - KS.SigmaS
    end
    KS.SigmaT =
        KS.SigmaS +
        sigmaA .* ones(size(KS.SigmaT)) +
        KS.convolutionmagnitude * ones(size(KS.SigmaT))

    # set estimate of contraction number
    L = KS.convolutionmagnitude / (1 / KS.dt + minimum(KS.SigmaT) + sqrt(2) / KS.dx)
    KS.contractionEstimate = L / (1 - L)

    # set moment order
    order = 0
    KS.DecompNumberElements = (order + 1) * (order + 1)

    # define total number of unknowns
    Nentries = KS.DecompNumberElements * KS.ny * KS.nx

    # decomposition order of kernel is Norder = (order+1)*(order+1)+1
    KS.DecompO, KS.DecompSigma, KS.DecompM = kerneldecomposition(KS, order) # Henyey Grenstein
    # precompute convolution matrix if needed
    if KS.convolutionmagnitude > 0
        KS.ScatteringKernelConvolution = computeScatteringKernelConvolution(KS)
    else
        KS.ScatteringKernelConvolution = zeros(KS.nquadpoints, KS.nquadpoints)
    end# precompute O*Sigma
    KS.OSigma = KS.DecompO * KS.DecompSigma

    # write Decomposition matrices for physical scattering
    M = KS.DecompM
    Sigma = KS.DecompSigma
    O = KS.DecompO
    OSigma = KS.OSigma

    # Set initial condition
    psi = deepcopy(KS.psi0)
    psi0 = deepcopy(KS.psi0)
    KS.lastpsi = deepcopy(KS.psi0)
    KS.psiC = deepcopy(KS.psi0)

    phiNew = zeros(KS.DecompNumberElements, KS.ny, KS.nx)
    @einsum phiNew[m, i, j] = M[m, q] * psi[q, i, j]

    rho = computeMoments(KS, psi)
    println("Initial mass = $(sum(rho[:])).")


    # define scattering kernel for source iteration
    S = deepcopy(KS.ScatteringKernelConvolution)
    for q = 1:size(S)[1]
        S[q, :] .*= KS.Q.weights[:]
    end

    KS.kernelTimesWeightsSparse = sparse(S)

    @showprogress 0.01 "Computing... " for k = 1:KS.nt

        # adapt dt if k*dt would be larger than tEnd
        dt = ifelse(k * KS.dt < KS.tEnd, KS.dt, KS.tEnd - k * KS.dt)

        # compute rhs
        rhsquads = sweeping(
            KS,
            deepcopy(KS.lastpsi),
            KS.lastpsi ./ KS.dt + KS.Source,
            KS.ScatteringKernelConvolution,
        )

        rhsmoms = zeros(KS.DecompNumberElements, KS.ny, KS.nx)
        @einsum rhsmoms[m, i, j] = M[m, q] * rhsquads[q, i, j]

        ## reshape to vector
        rrhsmoms = reshape(rhsmoms, Nentries)

        # define linear Map for Krylov solver
        D = LinearMap((x) -> linFunction(x, KS), Nentries, Nentries; ismutating = false) # muss in jeder Zeit iteration neu definiert werden, da KS sich aendert!

        # initial guess for Krylov solver
        #x0 = zeros(KS.DecompNumberElements,KS.nx,KS.ny)
        #@einsum x0[m,i,j] = M[m,q]*psi0[q,i,j]
        #x0 = reshape(x0,Nentries)

        # solve system with GMRES
        #println("############ starting krylov")
        println(sqrt(eps(real(eltype(rrhsmoms)))))
        phiNew = IterativeSolvers.gmres!(rrhsmoms, D, rrhsmoms, verbose = false)
        #println("done krylov")

        #phiNew = IterativeSolvers.gmres!(x0,D,rrhsmoms,verbose=false)
        phiNew = reshape(phiNew, (KS.DecompNumberElements, KS.nx, KS.ny))

        #visMoments(KS,phiNew,"OutputGMRES")

        # determine PsiNew
        OSigmaPhi = zeros(KS.nquadpoints, KS.nx, KS.ny)
        @einsum OSigmaPhi[q, i, j] = OSigma[q, moments] * phiNew[moments, i, j]

        # multiply with sigmaS
        @inbounds for j in KS.rangex, i in KS.rangey
            OSigmaPhi[:, i, j] = OSigmaPhi[:, i, j] .* KS.SigmaS[i, j]
        end

        # perform sweeping to update psi   
        psiNew = sweeping(
            KS,
            deepcopy(KS.lastpsi),
            KS.Source .+ KS.lastpsi / KS.dt + OSigmaPhi,
            KS.ScatteringKernelConvolution,
        )
        KS.lastpsi .= psiNew
    end

    rho = computeMoments(KS, KS.lastpsi)
    println("Final mass = $(sum(rho[:])).")

    vis(KS, KS.lastpsi)
    return rho
end

function solveWithSweeping!(KS::KineticSolver)

    # set max number of sweeping steps
    KS.nSweepingSteps = 1000
    KS.epsSweeping = 1e-3#1e-7;

    # update SigmaT if convolution is active
    if KS.ictype == "linesource"
        println("Problem Type is ", KS.ictype)
        sigmaA = 0.0 # hard coded for linesource
    else
        sigmaA = KS.SigmaT - KS.SigmaS
    end
    KS.SigmaT =
        KS.SigmaS +
        sigmaA .* ones(size(KS.SigmaT)) +
        KS.convolutionmagnitude * ones(size(KS.SigmaT))

    # set estimate of contraction number
    L =
        (KS.convolutionmagnitude + maximum(KS.SigmaS)) /
        (1 / KS.dt + minimum(KS.SigmaT) + sqrt(2) / KS.dx)
    KS.contractionEstimate = L / (1 - L)

    # compute artificial scattering kernel
    KS.ScatteringKernelConvolution = computeScatteringKernelConvolution(KS) # Convo

    # Set initial condition
    psi = zeros(KS.nquadpoints, KS.ny, KS.nx)
    psi = deepcopy(KS.psi0)
    KS.lastpsi = deepcopy(KS.psi0)
    rho = computeMoments(KS, psi)
    println("Initial mass = $(sum(rho[:])).")

    nicewrite(
        string(KS.outputfolderprefix, "data/source.txt"),
        computeMoments(KS, KS.Source),
    )

    # define scattering kernel for source iteration
    asKernel = deepcopy(KS.ScatteringKernelConvolution) # artificial inscattering
    physKernel = ones(size(asKernel)) ./ 4.0 / pi # add physical inscattering
    g = 0.01
    HG = computeScatteringKernelHenyeyGreenstein(KS, g)
    # decomposition order of kernel is Norder = (order+1)*(order+1)+1
    KS.DecompO, KS.DecompSigma, KS.DecompM = kerneldecomposition(KS, 1) # Henyey Grenstein


    KS.kernelTimesWeightsSparse = sparse(asKernel)

    #println("HG = ",HG)
    #println("OSigmaM = ",KS.DecompO*KS.DecompSigma*KS.DecompM)
    #println("wQ = ",KS.Q.weights[:])
    for q = 1:size(physKernel)[1]
        asKernel[q, :] .*= KS.Q.weights[:]
        physKernel[q, :] .*= KS.Q.weights[:]
    end
    physKernel = KS.DecompO * KS.DecompSigma * KS.DecompM

    @showprogress 0.01 "Computing... " for k = 1:KS.nt

        # adapt dt if k*dt would be larger than tEnd
        dt = ifelse(k * KS.dt < KS.tEnd, KS.dt, KS.tEnd - k * KS.dt)

        # first euler step
        psiNew = sweepingWithKernels(
            KS,
            deepcopy(KS.lastpsi),
            KS.Source .+ KS.lastpsi / KS.dt,
            asKernel,
            physKernel,
        )
        KS.lastpsi .= psiNew
        #vis(KS,psi)
        #return;
    end
    vis(KS, KS.lastpsi) # do all post computational plotting
    rho = computeMoments(KS, KS.lastpsi)
    println("End mass = $(sum(rho[:])).")
end # solve

function solveRSNWithSweeping!(KS::KineticSolver)

    # set max number of sweeping steps
    KS.nSweepingSteps = 1000
    KS.epsSweeping = 1e-3#1e-7;

    # update SigmaT if convolution is active
    if KS.ictype == "linesource"
        println("Problem Type is ", KS.ictype)
        sigmaA = 0.0 # hard coded for linesource
    else
        sigmaA = KS.SigmaT - KS.SigmaS
    end
    KS.SigmaT =
        KS.SigmaS +
        sigmaA .* ones(size(KS.SigmaT)) +
        KS.convolutionmagnitude * ones(size(KS.SigmaT))

    # set estimate of contraction number
    L =
        (KS.convolutionmagnitude + maximum(KS.SigmaS)) /
        (1 / KS.dt + minimum(KS.SigmaT) + sqrt(2) / KS.dx)
    KS.contractionEstimate = L / (1 - L)

    # compute artificial scattering kernel
    KS.ScatteringKernelConvolution = computeScatteringKernelConvolution(KS) # Convo

    # Set initial condition
    psi = zeros(KS.nquadpoints, KS.ny, KS.nx)
    psi = deepcopy(KS.psi0)
    KS.lastpsi = deepcopy(KS.psi0)
    rho = computeMoments(KS, psi)
    println("Initial mass = $(sum(rho[:])).")

    nicewrite(
        string(KS.outputfolderprefix, "data/source.txt"),
        computeMoments(KS, KS.Source),
    )

    # define scattering kernel for source iteration
    asKernel = deepcopy(KS.ScatteringKernelConvolution) # artificial inscattering
    physKernel = ones(size(asKernel)) ./ 4.0 / pi # add physical inscattering
    g = 0.01
    HG = computeScatteringKernelHenyeyGreenstein(KS, g)
    # decomposition order of kernel is Norder = (order+1)*(order+1)+1
    KS.DecompO, KS.DecompSigma, KS.DecompM = kerneldecomposition(KS, 1) # Henyey Grenstein


    KS.kernelTimesWeightsSparse = sparse(asKernel)

    #println("HG = ",HG)
    #println("OSigmaM = ",KS.DecompO*KS.DecompSigma*KS.DecompM)
    #println("wQ = ",KS.Q.weights[:])
    for q = 1:size(physKernel)[1]
        asKernel[q, :] .*= KS.Q.weights[:]
        physKernel[q, :] .*= KS.Q.weights[:]
    end
    physKernel = KS.DecompO * KS.DecompSigma * KS.DecompM

    @showprogress 0.01 "Computing... " for k = 1:KS.nt

        # adapt dt if k*dt would be larger than tEnd
        dt = ifelse(k * KS.dt < KS.tEnd, KS.dt, KS.tEnd - k * KS.dt)

        # first euler step
        psiNew = sweepingWithKernels(
            KS,
            deepcopy(KS.lastpsi),
            KS.Source .+ KS.lastpsi / KS.dt,
            asKernel,
            physKernel,
        )
        KS.lastpsi .= psiNew

        # perform rotation and interpolation step, or convolution step, or other modifications
        methods!(KS, KS.lastpsi, psiNew, k)
        #vis(KS,psi)
        #return;
    end
    vis(KS, KS.lastpsi) # do all post computational plotting
    rho = computeMoments(KS, KS.lastpsi)
    println("End mass = $(sum(rho[:])).")
end # solve

function solve!(KS::KineticSolver)
    #close("all")
    psi = zeros(KS.nquadpoints, KS.ny, KS.nx)
    psiLastTimeStep = zeros(KS.nquadpoints, KS.ny, KS.nx)
    flux = zeros(KS.nquadpoints, KS.ny, KS.nx)
    psi = deepcopy(KS.psi0) # Set initial condition
    rho = computeMoments(KS, psi)
    println("Initial mass = $(sum(rho[:])).")

    #nicewrite(string(KS.outputfolderprefix,"data/rho0.txt"),computeMoments(KS,psi))
    nicewrite(
        string(KS.outputfolderprefix, "data/source.txt"),
        computeMoments(KS, KS.Source),
    )

    summarizeAsPlot(KS, psi) #close("all")



    # precompute convolution matrix if needed
    if KS.convolutionmagnitude > 0
        KS.ScatteringKernelConvolution = computeScatteringKernelConvolution(KS)
    else
        KS.ScatteringKernelConvolution = zeros(KS.nquadpoints, KS.nquadpoints)
    end



    # low rank computation?
    if KS.rank > 0
        dA = zeros(KS.nquadpoints, (KS.ny) * (KS.nx))
        psireshaped = reshape(psi, KS.nquadpoints, KS.ny * KS.nx)
        s0 = svd(psireshaped)
        U0 = s0.U[1:end, 1:KS.rank]
        S0 = diagm(0 => s0.S[1:KS.rank])
        V0 = s0.V[1:end, 1:KS.rank]
        U1 = zeros(size(U0))
        V1 = zeros(size(V0))
        S1 = zeros(size(S0))
        K1 = zeros(size(U0 * S0))
        L1 = zeros(size(V0 * S1'))
        Y0 = zeros(KS.nquadpoints, (KS.ny) * (KS.nx))
    end

    @showprogress 0.01 "Computing... " for k = 1:KS.nt

        # force periodicity (if wanted)
        forcePeriodicity!(KS, psi)

        # adapt dt if k*dt would be larger than tEnd
        dt = ifelse(k * KS.dt < KS.tEnd, KS.dt, KS.tEnd - k * KS.dt)

        # store old psi to average later on
        #@inbounds for q=1:length(psi);	psiLastTimeStep[q] = psi[q]	end
        psiLastTimeStep .= psi

        # first euler step
        solveFlux!(KS, psi, flux)
        Euler!(KS, psi, flux, dt)

        # second euler step
        solveFlux!(KS, psi, flux)
        Euler!(KS, psi, flux, dt)

        #average
        psi .= (psi .+ psiLastTimeStep) ./ 2.0

        # add source contribution
        psi .+= KS.Source .* KS.dt

        # Note: So far, low rank and rSN DO NOT WORK TOGETHER. Either lowrank or rSN.
        if KS.rank > 0
            @inbounds for q = 1:length(psi)
                dA[q] = psi[q] - psiLastTimeStep[q]
            end
            K1 = U0 * S0 + dA * V0
            U1, S1 = qr(K1)
            U1 = Matrix(U1)
            S1 = S1 - U1' * dA * V0
            L1 = V0 * S1' + dA' * U1
            V1, S1 = qr(L1)
            V1 = Matrix(V1)
            S1 = S1'
            Y1 = U1 * S1 * V1'
            for i = 1:length(Y0), Y0[i] in Y1[i]
            end
            for i = 1:length(U0), U0[i] in U1[i]
            end
            for i = 1:length(V0), V0[i] in V1[i]
            end
            for i = 1:length(S0), S0[i] in S1[i]
            end
            for i = 1:length(psi), psi[i] in Y0[i]
            end
        else
            # perform rotation and interpolation step, or convolution step, or other modifications
            methods!(KS, psi, psiLastTimeStep, k)
        end
    end
    vis(KS, psi) # do all post computational plotting
end # solve
