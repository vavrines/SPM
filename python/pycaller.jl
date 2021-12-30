using PyCall

# workflow to import local python module
pushfirst!(PyVector(pyimport("sys")["path"]), @__DIR__)
pushfirst!(PyVector(pyimport("sys")["path"]), "$(@__DIR__)/PlatonicSolidsQuadratures")
asn = pyimport("asSN")

# ny, nx, cfl, points, weights, sigmaas, beta
ps = asn.initparams()

# asn.execute(asn.linesource(), ps)
