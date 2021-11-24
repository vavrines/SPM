using DelimitedFiles
using LinearAlgebra
using FastGaussQuadrature
include("Quadrature.jl")
include("standardQuadrature.jl")

# gauss
for n=2:2:100
	p,w = Quadrature(n,1)
	p = [p; p[:,1] p[:,2] -p[:,3]]
	w = [w;w]

	m = [p w]
	display(n)
	name = string(n)*"_gauss.txt"
	writedlm(name,m)
end
#
## octa
#for n=2:2:50
#	p,w = Quadrature(n,2)
#	m = [p w]
#	display(n)
#	name = string(n)*"_s_octa.txt"
#	writedlm(name,m)
#end
#
## ico
#for n=2:2:50
#	p,w = Quadrature(n,3)
#	m = [p w]
#	display(n)
#	name = string(n)*"_s_ico.txt"
#	writedlm(name,m)
#end
