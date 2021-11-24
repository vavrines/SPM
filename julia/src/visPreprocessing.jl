using PyPlot

function visPreprocessing(KS::KineticSolver, file::String)

    rho = readdlm(file, ',')
    #println("Final mass = $(sum(rho[:]))." ) 
    #nicewrite(string(KS.outputfolderprefix,"data/rhofinal.txt"),rho)


    if KS.ictype == "linesource"
        fig, ax = subplots(figsize = (10.5, 8), dpi = 100)
        minVal = 0.0
        maxVal = 0.5
        #	title("$(KS.Q.quadtype), nq=$(KS.nquadpoints), cs=$(KS.convolutionmagnitude), rot=$(KS.rotationmagnitude)", fontsize=20)
        pcolormesh(rho, cmap = "RdBu_r", vmin = minVal, vmax = maxVal)
        colorbar()
        PyPlot.savefig(string(KS.outputfolderprefix, "rho.png"))

        rhoana = readdlm("../vis/exactLineSource.txt", ',')
        #fig, ax = subplots(figsize=(10.5, 8), dpi=100)
        plot(rhoana[:, 1], rhoana[:, 2])
        rhong = rho[3:end-2, 3:end-2] # remove ghosts
        cuthori = rhong[Int64(KS.ny / 2), :]
        x = range(-1.5, stop = 1.5, length = KS.nx - 4)
        plot(x, cuthori)
        cutverti = rhong[:, Int64(KS.nx / 2)]
        plot(x, cutverti)
        cutdiag = [rhong[i, i] for i = 1:KS.nx-4]
        plot(sqrt(2) * x, cutdiag)
        cutdiag = [rhong[KS.ny-4+1-i, i] for i = 1:KS.nx-4]
        plot(sqrt(2) * x, cutdiag)
        title(
            "$(KS.Q.quadtype), nq=$(KS.nquadpoints), convolution=$(KS.convolutionmagnitude) ",
            fontsize = 30,
        )
        PyPlot.savefig(string(KS.outputfolderprefix, "rhocut.png"))

    elseif KS.ictype == "checkerboard"
        fig, ax = subplots(figsize = (10.5, 8), dpi = 100)
        minVal = -7.0
        maxVal = 0.0
        title(
            "$(KS.Q.quadtype), nq=$(KS.nquadpoints), cs=$(KS.convolutionmagnitude)",
            fontsize = 20,
        )
        pcolormesh(log10.(abs.(rho)), cmap = "RdBu_r", vmin = minVal, vmax = maxVal)
        colorbar()
        PyPlot.savefig(string(KS.outputfolderprefix, "rho_log10.png"))

    elseif KS.ictype == "beam"

        fig, ax = subplots(figsize = (10.5, 8), dpi = 100)
        minVal = -7.0
        maxVal = 0.0
        title(
            "$(KS.Q.quadtype), nq=$(KS.nquadpoints), cs=$(KS.convolutionmagnitude)",
            fontsize = 20,
        )
        pcolormesh(log10.(abs.(rho[3:end-2, 3:end-2])), cmap = "RdBu_r")#,vmin=minVal, vmax=maxVal)
        colorbar()
        PyPlot.savefig(string(KS.outputfolderprefix, "rho_log10.png"))

        fig, ax = subplots(figsize = (10.5, 8), dpi = 100)
        minVal = -7.0
        maxVal = 0.0
        title(
            "$(KS.Q.quadtype), nq=$(KS.nquadpoints), cs=$(KS.convolutionmagnitude)",
            fontsize = 20,
        )
        pcolormesh(((rho[3:end-2, 3:end-2])), cmap = "RdBu_r")#,vmin=minVal, vmax=maxVal)
        colorbar()
        PyPlot.savefig(string(KS.outputfolderprefix, "rho.png"))


    end
end
