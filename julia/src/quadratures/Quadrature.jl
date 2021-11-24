include("geometryhelper.jl")
include("icosahedronQuadrature.jl")
include("octaQuadrature.jl")
include("standardQuadrature.jl")
include("dim2Quadrature.jl")
using FastGaussQuadrature


export Quadrature
export rotateandinterpolate!
export rotateOnly!

mutable struct Quadrature

	norder::Int64
	nquadpoints::Int64
	pointsxyz::Array{Float64,2}
	pointsmuphi::Array{Float64,2}
    weights::Array{Float64,1}

    muphiOLD::Array{Float64,2}
    xyzOLD::Array{Float64,2}

    triangulation::Array{Int64,2}

    interpolationids::Array{Int64,2}
	interpolationweights::Array{Float64,2}
	neighbors::Array{Int64,2}

    rotationmatrix::Array{Float64,2}

	quadtype::String;

	function Quadrature()
		return Quadrature(5,1)
	end

    function Quadrature(N::Int64,quadtype=1::Int)
		n = N
		quadTypeString = "quad"

		if quadtype==1 # standard quad
			xyz,weights = createPointStandards(n)
			triangulation = zeros(Int,1,1)
			neighbors = zeros(Int64,length(weights),6)
			nquadpoints = size(xyz)[1]	
			quadTypeString = "Tensorized"		

		elseif quadtype == 2 # wn octa 
			xyz,triangulation =  createPointsOcta(n)
			nquadpoints = size(xyz)[1]			
			weights = createWeights(nquadpoints,xyz,triangulation)
			neighbors = zeros(Int64,length(weights),6)
			for i=1:length(weights)
				tmp = (1 .+ rem.(findall(triangulation[:].==i) .- 1,size(triangulation,1)))
				neighbors[i,1:length(tmp)] = tmp;
			end
			quadTypeString = "Octaeder"	
		elseif quadtype == 3 # wn ico
			xyz,triangulation =  createPointsIco(n)	
			nquadpoints = size(xyz)[1]			
			weights = createWeights(nquadpoints,xyz,triangulation)
			neighbors = zeros(Int64,length(weights),6)
			for i=1:length(weights)
				tmp = (1 .+ rem.(findall(triangulation[:] .==i) .- 1,size(triangulation,1)))
				neighbors[i,1:length(tmp)] = tmp;
			end
			quadTypeString = "Icosaeder"	
		elseif quadtype == 4 # 2d
			xyz,weights = createPoint2d(n)
			triangulation = zeros(Int,1,1)
			neighbors = zeros(Int64,length(weights),6)
			nquadpoints = size(xyz)[1]
			quadTypeString = "test_2d"
		elseif quadtype == 5 # from file
			println("Please enter the quadrature file to be loaded and press enter.")
			println("(empty field + enter = quadrature.txt)")
			file = readline()
			if isempty(file)
				file = "quadrature.txt"
			end
			println("Trying to read $file.")
			M = readdlm(file,',')
			xyz = M[:,1:3];
			weights = M[:,4];
			triangulation = zeros(Int,1,1)
			neighbors = zeros(Int64,length(weights),6)
			nquadpoints = size(xyz)[1]
			quadTypeString = "fromfile"
		end
		
		muphi = zeros(Float64,nquadpoints,2)
		xyz2muphi!(xyz,muphi)

		interpolationids = zeros(Int64,nquadpoints,3)
		interpolationweights = zeros(Float64,nquadpoints,3)

		rotationmatrix = zeros(3,3); 
		rotationmatrix[1,1] = 1.0
		rotationmatrix[2,2] = 1.0
		rotationmatrix[3,3] = 1.0

		xyzOLD = deepcopy(xyz)
		muphiOLD = deepcopy(muphi)


		new(n,nquadpoints,xyz,muphi,weights,muphiOLD,xyzOLD,triangulation,
			interpolationids,interpolationweights,neighbors,rotationmatrix,quadTypeString)
    end
end # struct

function createWeights(n::Int64,xyz::Array{Float64,2},triangles::Array{Int64,2})
	weights = zeros(n)
	#Compute the quadrature weights for each node
	#Therefore we iterate over all triangles and increase the weight
	#of each of the nodes of the triangle by the area of the triangle that
	#is assigned to it
	

	weights = zeros(n)
	nTriangles = size(triangles)[1]
	xy = zeros(Float64,3)
	yz = zeros(Float64,3)
	zx = zeros(Float64,3)
	mid = zeros(Float64,3)
	
	for n=1:nTriangles

		# get three nodes of a triangle
		i,j,k = triangles[n,:]

		# get the corresponding points 
		x = xyz[i,:]
		y = xyz[j,:]
		z = xyz[k,:]

		# Now get the midpoint of the triangle and the midpoints along the lines
		mid = (x+y+z)/3.;
		xy = (x+y)/2.;
		yz = (y+z)/2.;
		zx = (z+x)/2.;
 
		# These points still have to be projected onto the sphere
		mid = mid/norm(mid,2);
		xy = xy/norm(xy,2);
		yz = yz/norm(yz,2);
		zx = zx/norm(zx,2);

		#By these four points, plus x,y,z we can span 6 triangles
		#Each of these triangles is assigned to one of the three nodes of the triangle x,y,z.
		#the area = the weight
		
		# for i
		weights[i] += getarea(x,mid,xy,"sphere")
		weights[i] += getarea(x,mid,zx,"sphere")

		# for j
		weights[j] += getarea(y,mid,xy,"sphere")
		weights[j] += getarea(y,mid,yz,"sphere")


		# for k
		weights[k] += getarea(z,mid,yz,"sphere")
		weights[k] += getarea(z,mid,zx,"sphere")


	end

    return weights
end

function applyinterpolation!(phiNew::Array{Float64,3},phi::Array{Float64,3},Q::Quadrature)
	nq,ni,nj = size(phi)
	#@parallel
		
	for j=1:nj
			for i=1:ni
				for q=1:nq
					id1 = Q.interpolationids[q,1]
					id2 = Q.interpolationids[q,2]
					id3 = Q.interpolationids[q,3]
			
					w1 = Q.interpolationweights[q,1]
					w2 = Q.interpolationweights[q,2]
					w3 = Q.interpolationweights[q,3]
					
					phiNew[q,i,j] = w1*phi[id1,i,j]+w2*phi[id2,i,j]+w3*phi[id3,i,j]
				
			end
		end
	end
	


	# fix to be conservative
	for j=1:nj
		for i=1:ni
			newnorm = 0.0;
			oldnorm = 0.0
			for q=1:nq
				newnorm += phiNew[q,i,j]*Q.weights[q]
				oldnorm += phi[q,i,j]*Q.weights[q]
			end
			if oldnorm>0.0
				for q=1:nq
					phiNew[q,i,j] *= (oldnorm/newnorm)
				end
			end
		end
	end
end




#Get the ID of the triangle to which v belongs
#if v belongs to multiple triangles because it is on an edge/vertix this
#returns still just one id as the interpolation would be the same,
#independent of the triangle id returned.
#If no triangle could be assigned, we throw an error
function gettriangleid(v::Array{Float64,1},Q::Quadrature,neighbors::Array{Int64,1})
	# this gives the rows of the triangulation matrix that contain i
	# check only neighbors first, usually this is enough
	for i=1:length(neighbors)
		if neighbors[i]>0 # the poles only have 4, not 6 neighbors
			a = Q.xyzOLD[Q.triangulation[neighbors[i],1],:]
			b = Q.xyzOLD[Q.triangulation[neighbors[i],2],:]
			c = Q.xyzOLD[Q.triangulation[neighbors[i],3],:]
		
			if isprojectedpointinsidetriangle(v,a,b,c)
				#We still have got to make sure that we don't pick the triangle from the opposite side.
				#Therefore it should be sufficient to check what the distance of the point
				#to one of the vertices of the triangle is.
				if norm(a-v)<1.0
					return neighbors[i];
				end
			end
		end
	end
	
	# check everything
	for i=1:size(Q.triangulation,1)
		a = Q.xyzOLD[Q.triangulation[i,1],:]
		b = Q.xyzOLD[Q.triangulation[i,2],:]
		c = Q.xyzOLD[Q.triangulation[i,3],:]
	
		if isprojectedpointinsidetriangle(v,a,b,c)
			#We still have got to make sure that we don't pick the triangle from the opposite side.
			#Therefore it should be sufficient to check what the distance of the point
			#to one of the vertices of the triangle is.
			if norm(a-v)<1.0
				return i;
			end
		end
	end
	println(v)
	println(norm(v)-1.0)
	error("We should not reach this part!")
end



function prepareinterpolation!(Q::Quadrature)
	# find a traingle id
	for i=1:Q.nquadpoints
		triangleID = gettriangleid(Q.pointsxyz[i,:],Q,Q.neighbors[i,:])
		Q.interpolationids[i,:] = Q.triangulation[triangleID,:]
		v = Q.pointsxyz[i,:]
		a = Q.xyzOLD[Q.interpolationids[i,1],:]
		b = Q.xyzOLD[Q.interpolationids[i,2],:]
		c = Q.xyzOLD[Q.interpolationids[i,3],:]
		
		Q.interpolationweights[i,:] = computeinterpolationweights(v,a,b,c)';
	end
end


function rotate!(Q::Quadrature)
	for i=1:size(Q.pointsxyz,1)
		Q.pointsxyz[i,:] = (Q.rotationmatrix*(Q.pointsxyz[i,:]))';
		Q.pointsxyz[i,:] = Q.pointsxyz[i,:]/norm(Q.pointsxyz[i,:])
	end
	xyz2muphi!(Q.pointsxyz,Q.pointsmuphi);

end


	

function rotateandinterpolate!(alpha::Float64,beta::Float64,gamma::Float64,
	Q::Quadrature,phiNew::Array{Float64,3},phi::Array{Float64,3},direction::String)

	getrotationmatrix!(Q.rotationmatrix,alpha,beta,gamma,direction)
	rotate!(Q)
	prepareinterpolation!(Q)
	applyinterpolation!(phiNew,phi,Q)
end

function rotateOnly!(u::Array{Float64,1},theta::Float64,
	Q::Quadrature,direction::String)
	getrotationmatrix!(Q.rotationmatrix,u,theta,direction)
	rotate!(Q)
	Q.xyzOLD = deepcopy(Q.pointsxyz)
	Q.muphiOLD = deepcopy(Q.pointsmuphi)

end


function rotateandinterpolate!(u::Array{Float64,1},theta::Float64,
	Q::Quadrature,phiNew::Array{Float64,3},phi::Array{Float64,3},direction::String)
	getrotationmatrix!(Q.rotationmatrix,u,theta,direction)
	rotate!(Q)
	prepareinterpolation!(Q)
	applyinterpolation!(phiNew,phi,Q)
	Q.xyzOLD = deepcopy(Q.pointsxyz)
	Q.muphiOLD = deepcopy(Q.pointsmuphi)

end
