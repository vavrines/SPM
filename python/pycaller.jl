using PyCall
np = pyimport("numpy")

# workflow to import local python module
pushfirst!(PyVector(pyimport("sys")["path"]), @__DIR__)
pushfirst!(PyVector(pyimport("sys")["path"]), "$(@__DIR__)/PlatonicSolidsQuadratures")
asn = pyimport("asSN")
pyimport("importlib")["reload"](asn) # hot update

# ny, nx, cfl, points, weights, sigmaas, beta
ps = asn.initparams(50, 50, 5)
# y0, y1, x0, x1, t, fsigmas, fsigmaa, fpsi0, fsource, style
p2s = asn.linesource(); p2s[5] = 0.1 # t
sols = asn.initsol(p2s, ps)
psi0 = deepcopy(sols[1])

#psi0, dt, dx, dy, nt, qpoints, qweights, sigmas, sigmaa, sigmaas, beta, source
psi = asn.run(sols...); sols[1] .= psi

# plot
begin
    using Plots
    
    y0, y1, x0, x1 = p2s[1:4]
    dx, dy = sols[3:4]
    qweights = sols[7]

    ρ = zeros(ps[1], ps[2])
    for i in axes(ρ, 1), j in axes(ρ, 2)
        ρ[i, j] = sum(qweights .* psi[:, i+2, j+2])
    end

    contourf(x0:dx:x1, y0:dy:y1, ρ)
end
