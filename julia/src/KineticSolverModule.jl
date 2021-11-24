module KineticSolverModule

import IterativeSolvers

include("MyIOModule.jl")
using .MyIOModule

using ProgressMeter, DelimitedFiles, LinearAlgebra, Dates, Revise, PyPlot
using Einsum, Arpack, LinearMaps, SparseArrays
using DelimitedFiles, SphericalHarmonics, TypedPolynomials, GSL

export KineticSolver, computeMoments, corysolve!, solve!, sweeping
export setTestCaseInitialCondition!,
    computeHenyeyGreensteinScatteringMatrix,
    solveWithSweeping!,
    solveRSNWithSweeping!,
    sweepingWithKernels,
    vis,
    methods!
export kerneldecomposition, metadecomposition, kerneldecomposition
export computeScatteringKernelHenyeyGreenstein, computeScatteringKernelConvolution
export checkexpansionHG, solveWithSweeping!, computeScatteringKernelConvolution

include("quadratures/Quadrature.jl")
include("problems/getProblemSpecificParameters.jl")

mutable struct KineticSolver
    # domain is given by [x0,x1] x [y0,y1]
    x0::Float64
    x1::Float64
    y0::Float64
    y1::Float64
    dx::Float64 # corresponding equidistant grid spacing
    dy::Float64 # same for y
    nx::Int64 # number of cells +2 ghost cells
    ny::Int64 # same for y

    norder::Int64 # Order of the quadrature
    nquadpoints::Int64 # Resulting number of quadraturepoints
    #for SN: nquadpoints = norder^2, but for other quadratures this is not true

    # index array for the specified variables
    # later we just want to write for i=rangex instead of for i=1:N
    rangex::Array{Int64,1}
    rangey::Array{Int64,1}
    rangequad::Array{Int64,1}

    ThisQuadratureType::String # Name of the Quadrature being used
    Q::Quadrature              # Quadrature struct

    rotationmagnitude::Float64 # if rotation is chosen
    convolutionmagnitude::Float64 # convolution strength
    convolutionwidth::Float64 # controls how isotropically scattering happens
    convolutiontoscattering::Float64 # controls convolution strength in regions with high scattering

    cfl::Float64 # cfl number
    dt::Float64 # time step size
    tEnd::Float64  # time until we simulate
    nt::Int64 # number of time steps

    tmp::Array{Float64,1}
    periodicXflag::Bool # periodic x boundary
    periodicYflag::Bool
    dim1flag::Bool # perform 1d simulation
    rank::Int64

    SigmaS::Array{Float64,2} # Matrix of Scattering crossection
    SigmaT::Array{Float64,2} # Matrix of Total crossection

    ictype::String ## Name of Initial Condition
    Source::Array{Float64,3} # Tensor for time independent source
    psi0::Array{Float64,3} # Tensor for the phi(t=0,x,y,q)

    outputfolderprefix::String # folder name prefix for storing the data
    staticCounter::Int

    nSweepingSteps::Int
    epsSweeping::Float64

    DecompO::Array{Float64,2}
    DecompSigma::Array{Float64,2}
    DecompM::Array{Float64,2}
    OSigma::Array{Float64,2}
    lastpsi::Array{Float64,3}
    DecompNumberElements::Int64

    ScatteringKernelConvolution::Array{Float64,2}
    OutscatteringVector::Array{Float64,1}

    contractionEstimate::Float64
    kernelTimesWeightsSparse::SparseMatrixCSC{Float64,Int64}
    psiC::Array{Float64,3} # Tensor for the phi(t=0,x,y,q) for one shot 

    function KineticSolver(configfilename::String)
        # read the following words from text file
        allowed = [
            "nx",
            "ny",
            "cfl",
            "quadraturetype",
            "quadratureorder",
            "convolutionmagnitude",
            "convolutionwidth",
            "convolutiontoscattering",
            "testcaseid",
            "rotationmagnitude",
            "periodicXflag",
            "periodicYflag",
            "dim1flag",
            "whichrank",
        ]
        D = readtodict(configfilename, allowed)

        #				 read data from config file start					  #
        nx = Int64(D["nx"])
        ny = Int64(D["ny"])
        norder = Int64(D["quadratureorder"])
        Q = Quadrature(norder, Int64(D["quadraturetype"]))
        rotationmagnitude = D["rotationmagnitude"]
        rank = Int64(D["whichrank"])
        cfl = D["cfl"]
        PossibleQuadratureTypes = ["tensorized", "octa", "ico", "dim2", "fromfile"]
        PossibleICTypes = ["linesource", "checkerboard", "dim1", "beam"]
        ThisQuadratureType = PossibleQuadratureTypes[Int64(D["quadraturetype"])]
        periodicXflag = Int64(D["periodicXflag"]) == 1 ? true : false
        periodicYflag = Int64(D["periodicYflag"]) == 1 ? true : false
        dim1flag = Int64(D["dim1flag"]) == 1 ? true : false
        ThisICType = PossibleICTypes[Int64(D["testcaseid"])]
        if dim1flag == true ## 1D FLAG FORCES THESE THINGS:
            ny = 3
            periodicYflag = true
            ThisICType = "dim1"
            text = "You chose a 1d simulation, this forces: ny = 1, periodicYflag = true, ThisICType = 3, ThisQuadratureType = 4"
            printstyled(text, color = :red)
        end
        #some variables can be deduced immediately though not being part of the config
        nquadpoints = Q.nquadpoints
        rangex = collect(1:nx) .+ 2 #range for x and y considering ghost cells
        rangey = collect(1:ny) .+ 2
        rangequad = collect(1:nquadpoints)
        convolutionmagnitude = Float64(D["convolutionmagnitude"])
        convolutionwidth = Float64(D["convolutionwidth"])
        convolutiontoscattering = Float64(D["convolutiontoscattering"])
        #				 read data from config file stop					  #

        ######## Get problem specific parameters start  ########
        x0, x1, y0, y1, tEnd, SigmaS, SigmaT, Source, psi0 =
            getProblemSpecificParameters(nx, ny, nquadpoints, ThisICType, Q.pointsxyz)
        # and deduce other variables from that
        dx = (x1 - x0) / nx
        dy = (y1 - y0) / ny
        dt = cfl / 2 * (dx * dy) / (dx + dy)
        nt = ceil(tEnd / dt)
        tmp = zeros(Float64, nquadpoints) # tmp array for in place mult
        ######## Get problem specific parameters stop  ########

        identifier = string(Dates.now(), "/")
        outputfolderprefix = string("../out/", identifier)
        mkdir(outputfolderprefix)
        mkdir(string(outputfolderprefix, "data/"))
        cp("config.txt", string(outputfolderprefix, "config.txt"))
        nicewrite(string(outputfolderprefix, "data/sigmaS.txt"), SigmaS)
        nicewrite(string(outputfolderprefix, "data/sigmaT.txt"), SigmaT)
        #nicewrite(string(outputfolderprefix,"data/sourceAngular.txt"),Source)
        #nicewrite(string(outputfolderprefix,"data/phi0.txt"),phi0)
        nicewrite(string(outputfolderprefix, "data/quadpoints.txt"), Q.pointsxyz)
        nicewrite(string(outputfolderprefix, "data/quadweights.txt"), Q.weights)


        # incomplete initialization:  weightMatrix will be computed at run time
        new(
            x0,
            x1,
            y0,
            y1,
            dx,
            dy,
            nx + 4,
            ny + 4,# spatial
            norder,
            nquadpoints,             # quadrature
            rangex,
            rangey,
            rangequad,        # range
            ThisQuadratureType,
            Q,          # quadrature
            rotationmagnitude,              # how much rotation
            convolutionmagnitude,           # how much convolution
            convolutionwidth,
            convolutiontoscattering,        # how much convolution in scattering/absorption regions
            cfl,
            dt,
            tEnd,
            nt,                 # time
            tmp,                            # tmp array for in place mult
            periodicXflag,
            periodicYflag, # are boundaries periodic?
            dim1flag,                     # 1d
            rank,                         # when performing low rank computations
            SigmaS,
            SigmaT,                # Sigma
            ThisICType,
            Source,
            psi0,       # IC/Q/PHI0
            outputfolderprefix,           # folder
            0,
        )
    end # constructor
end # struct

include("vis.jl")
include("methods.jl")
include("sweeping.jl")
include("kerneldecomposition.jl")
include("solves.jl")
include("computeKernels.jl")


end # module
