function createPoint2d(n::Int64)
    points,weights  = computeXYZandWeightsdim2(n)
	return points,weights
end


function computeXYZandWeightsdim2(_norder::Int64)
        n = _norder
        pointsxyz = zeros(n,3); # Even though we only need the x and y coordinate 
        weights = zeros(n);
    
            
        # around z axis equidistant
		phi = [k/n*2*pi for k=0:n-1]

    
            
		for i=1:n
			pointsxyz[i,1] = cos(phi[i]) 
			pointsxyz[i,2] = sin(phi[i])
			pointsxyz[i,3] = 0.0
            
			weights[i] = 4*pi/n
		end
        
        return pointsxyz, weights
    end
