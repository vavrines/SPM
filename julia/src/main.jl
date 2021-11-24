cd(@__DIR__)
include("MyIOModule.jl")
include("KineticSolverModule.jl")

using .KineticSolverModule
using .MyIOModule
using DelimitedFiles

function main() # Executes all possible simulations resulting from the specified parameters below

    close("all")

    NX = [200]#[200]        # number of spatial cells in x and y direction
    CFL = [2.0]#[2.0 5.0 10 15 20]
    ORDER = [2, 3, 4, 5, 6, 7, 8, 9, 10]        # order of quadrature set
    ROTATION = [0]        # select rotation for rSN, 0 rotation deactivates rSN and uses SN
    QTYPE = [3]        # Type must be 1 for "tensorized" or 2 for  "octa" and 3 for "ico".
    CONVOLUTIONMGNTD = [0 7] #[0 0.5 1 1.5 2 3 4 5 6 7 8 9 10 12 14 16 18 20]        # select magnitude of convolution, 0 turns feature off
    CONVOLUTIONWIDTH = [4] #[2 2.5 3 3.5 4 4.5 5 6 7 8 9 10 12 14 16 18 20]        # select width of the convolution kernel, inf = isotropic scattering
    CONVOLUTIONTOSCAT = [0]        # setting this to anything but zero reduces conv. in dep. of total c.s.
    PROBLEM = [1]        # Linesource = 1; Checkerboard = 2
    PERIODICXFLAG = [0]        # shall the x domain be periodic?
    PERIODICYFLAG = [0]        # shall the y domain be periodic?
    DIM1FLAG = [0]        # shall the computation be a one dimensional (only x) computation?
    RANK = [0]        # select rank for low rank computation, 0 turns feature off

    counter = 1
    runtimes = zeros(length(CONVOLUTIONMGNTD) * length(ORDER))
    infoOrder = zeros(length(CONVOLUTIONMGNTD) * length(ORDER))
    infoMagnitude = zeros(length(CONVOLUTIONMGNTD) * length(ORDER))
    for nx in NX,
        cfl in CFL,
        order in ORDER,
        rotation in ROTATION,
        qtype in QTYPE,
        problem in PROBLEM,
        periodicXflag in PERIODICXFLAG,
        periodicYflag in PERIODICYFLAG,
        dim1flag in DIM1FLAG,
        rank in RANK,
        convMagn in CONVOLUTIONMGNTD,
        convWidth in CONVOLUTIONWIDTH,
        cts in CONVOLUTIONTOSCAT

        matchwords = [
            "nx",
            "ny",
            "cfl =",
            "quadratureorder =",
            "rotationmagnitude =",
            "quadraturetype = ",
            "convolutionmagnitude",
            "testcaseid",
            "periodicXflag",
            "periodicYflag",
            "dim1flag",
            "whichrank",
            "convolutionwidth",
            "convolutiontoscattering",
        ]

        replaceby = [
            "nx = $nx",
            "ny = $nx",
            "cfl = $cfl",
            "quadratureorder = $order",
            "rotationmagnitude = $rotation",
            "quadraturetype = $qtype",
            "convolutionmagnitude = $convMagn",
            "testcaseid = $problem",
            "periodicXflag = $periodicXflag",
            "periodicYflag = $periodicYflag",
            "dim1flag = $dim1flag",
            "whichrank = $rank",
            "convolutionwidth = $convWidth",
            "convolutiontoscattering = $cts",
        ]

        replacelineswithstring("config.txt", "config.txt", matchwords, replaceby)

        # Run
        KS = KineticSolver("config.txt")
        #@time solveWithSweeping!(KS)
        #@time solve!(KS)
        #@time solveRSNWithSweeping!(KS)
        start = time()
        @time corysolve!(KS)
        runtimes[counter] = time() - start
        infoOrder[counter] = order
        infoMagnitude[counter] = convMagn
        println("Run $counter finished.")
        counter = counter + 1

        close("all")
    end
    writedlm("runtimes", runtimes)
    writedlm("infoOrder", infoOrder)
    writedlm("infoMagnitude", infoMagnitude)
end

main()
