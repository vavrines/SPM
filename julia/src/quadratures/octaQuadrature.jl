function createPointsOcta(n::Int64)
    points,triangulation  = computexyz(n)
	return points,triangulation
end

function linspace(start,stop,n)
	return collect(range(start,stop=stop,length=n))
end

function computexyz(n::Int64,slerpflag=true::Bool)
	pt0 = [0.0, 0.0, 1.0]
	pt1 = [0.0, 1.0, 0.0]
	pt2 = [1.0, 0.0, 0.0]


	if slerpflag
		# use slerp instead of linspace
		pts01 = slerp(pt0,pt1,n)
		pts02 = slerp(pt0,pt2,n)
	else
		pts01 = linspace(pt0,pt1,n) #type is: n-element LinSpace{Array{Float64,1}}:
		pts02 = linspace(pt0,pt2,n)
	end



	
	nptsoctant = Int64(n*(n+1)/2)
	pts = zeros(3,nptsoctant)

	# Creates points in planar geometry
	counter = 0;
	for i=1:n
		if slerpflag
			if i==1
				tmp = pts01[:,1]
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



	# Now we have the quadrature points and triangles for a single octant
	ptsAll = deepcopy(pts);
	#mat ptsAll = pts; // the quadrature points for all octants, but with duplicated #entries

	tmp = deepcopy(pts); tmp[1,:] *= -1.0
	ptsAll = deepcopy(hcat(ptsAll,tmp));

	tmp = deepcopy(pts); tmp[2,:] *= -1.0;
	ptsAll = deepcopy(hcat(ptsAll,tmp));

	tmp = deepcopy(pts); tmp[1,:] *= -1.0; tmp[2,:] *= -1.0;
	ptsAll = deepcopy(hcat(ptsAll,tmp));

	tmp = deepcopy(ptsAll); tmp[3,:] *= -1.0;
	ptsAll = deepcopy(hcat(ptsAll,tmp));


	trianglesAll = deepcopy(triangles);
	for i=1:7
		trianglesAll = (hcat(trianglesAll,triangles .+i .*nptsoctant))
	end
	ptsAll,triangulation =unique(ptsAll,trianglesAll)

	ptsAll = permutedims(ptsAll)
	triangulation = permutedims(triangulation)
	return ptsAll, triangulation
end





