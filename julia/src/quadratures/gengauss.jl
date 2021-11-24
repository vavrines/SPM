using DelimitedFiles
using FastGaussQuadrature
include("standardQuadrature.jl")


for n=2:2:200
	p,w = createPointStandards(n)
	m = [p w]
	display(n)
	name = string(n)*"gauss.txt"
	writedlm(name,m)
end
