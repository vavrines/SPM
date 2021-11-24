
function createPointsIco(n::Int64)
    points,triangulation  = computeXYZico(n)
	return points,triangulation
end
function computeXYZico(n::Int64)
    # the vertices and faces of the regular icosahedron
    # with vertices on the unit sphere are hard coded
    # and stored in txt files
    vertices = readdlm("quadratures/icosahedron_vertices.txt",',',Float64)
    faces = readdlm("quadratures/icosahedron_faces.txt",',',Int)

    nfaces = 20;
    nvertices = 12;
    #close("all")

    n1,n2,n3 = faces[1,1],faces[1,2],faces[1,3]
    v1,v2,v3 = vertices[n1,:],vertices[n2,:],vertices[n3,:]

    p,t = getNodesAndTriangulationForSingleTriangle(n,v1,v2,v3)
   # plot3D(p[1,:],p[2,:],p[3,:],"x")

    ptsAll = deepcopy(p);
    trianglesAll = deepcopy(t);



    for f=2:20
        n1,n2,n3 = faces[f,1],faces[f,2],faces[f,3]
        v1,v2,v3 = vertices[n1,:],vertices[n2,:],vertices[n3,:]

        p,t = getNodesAndTriangulationForSingleTriangle(n,v1,v2,v3)
        nperface = size(p,2)
    #    plot3D(p[1,:],p[2,:],p[3,:],"x")
        ptsAll = deepcopy(hcat(ptsAll,p))
        trianglesAll = deepcopy(hcat(trianglesAll,t.+(f-1)*nperface))
    end
    #println(size(ptsAll))
    #println(size(trianglesAll))

    ptsAll,trianglesAll =unique(ptsAll,trianglesAll)



    # there will be n*(n+1)/2*20-30*n+1*12 points
	# and n*(n-2)*20+20 triangles
	ptsAll = permutedims(ptsAll)
	trianglesAll = permutedims(trianglesAll)
	return ptsAll,trianglesAll


end


function getNodesAndTriangulationForSingleTriangle(n::Int64,
    pt0::Array{Float64,1},pt1::Array{Float64,1},pt2::Array{Float64,1},slerpflag=true::Bool)


	if slerpflag
		# use slerp instead of linspace
		pts01 = slerp(pt0,pt1,n,true)
		pts02 = slerp(pt0,pt2,n,true)
	else
		pts01 = linspace(pt0,pt1,n) #type is: n-element LinSpace{Array{Float64,1}}:
		pts02 = linspace(pt0,pt2,n)
	end


	nptsoctant = Int64(n*(n+1)/2)
	pts = zeros(Float64,3,nptsoctant)

	# Creates points in planar geometry
	counter = 0;
	for i=1:n
		if slerpflag
			if i==1
				tmp = pts01[:,1]
			elseif i==n
				tmp = slerp(pts01[:,i],pts02[:,i],i,true)
			else
				tmp = slerp(pts01[:,i],pts02[:,i],i)
			end
		else
			tmp =linspace(pts01[i],pts02[i],i)
		end

		for j=1:i
			counter += 1
			if slerpflag
				pts[:,counter] = tmp[:,j]
			else
				pts[:,counter] = tmp[j]
			end
		end

	end
	# Project points onto sphere
	for i=1:nptsoctant
    	pts[:,i] = pts[:,i]/norm(pts[:,i])
	end

	# We now want to get the connectivity. Therefore enumerate all
	#the points (overkill from a computational point of view) and then
	#write the correct connectivity.

    # Matrix that assigns an ID to the points
   ids=zeros(Int64,n,n); # Too large, but not important.
   #Matrix that will later contain all triangles
   nTrianglesOctant = Int64(n*n-2*n+1);
   triangles = zeros(Int64,3,nTrianglesOctant);

	counter = 0
	for i=1:n
		for j=1:i
			counter +=1
			ids[i,j] = counter;
		end
	end

	# Now create triangles
	counter = 0
	tmp = zeros(Int64,1)
	for i=1:n
		for j=1:i-1
			tmp  = [ids[i,j],ids[i,j+1],ids[i-1,j]]
			counter += 1
			triangles[:,counter] = tmp;
		end
		if i<n
			for j=1:i-1
				tmp = [ids[i,j],ids[i,j+1],ids[i+1,j+1]]
				counter +=1
				triangles[:,counter] = tmp;
			end
		end
    end
    return pts,triangles;

end
