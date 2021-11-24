function getProblemSpecificParameters(nx, ny, nquadpoints, ThisICType, pointsxyz)

    if ThisICType == "linesource"
        SigmaS = 1 * ones(ny + 4, nx + 4)
        SigmaA = 0 * ones(ny + 4, nx + 4)
        SigmaT = SigmaS + SigmaA
        x0 = -1.5
        x1 = 1.5
        y0 = -1.5
        y1 = 1.5
        dx = (x1 - x0) / (nx)
        dy = (y1 - y0) / (ny)
        tEnd = 1.0
        phi0 = setInitialConditionsGaussian(nx, ny, nquadpoints, dx, dy, x0, y0)
        Source = zeros(Float64, nquadpoints, ny + 4, nx + 4)
    elseif ThisICType == "beam"
        SigmaS = 0.0 * ones(ny + 4, nx + 4)
        SigmaA = 0.0 * ones(ny + 4, nx + 4)
        SigmaT = SigmaS + SigmaA
        x0 = -1.5
        x1 = 1.5
        y0 = -1.5
        y1 = 1.5
        dx = (x1 - x0) / (nx)
        dy = (y1 - y0) / (ny)
        tEnd = 10.0

        beamdirection = [1 / sqrt(2), 1 / sqrt(2), 0]
        beamdirection = [0, 1, 0]
        Source = setInitialConditionsBeam(
            nx,
            ny,
            nquadpoints,
            dx,
            dy,
            x0,
            y0,
            pointsxyz,
            beamdirection,
        )
        phi0 = zeros(Float64, nquadpoints, ny + 4, nx + 4)
        v = 10.0
        for j = 3:nx+2
            for i = 3:ny+2
                y = y0 + dy / 2 + (i - 3) * dy
                x = x0 + dx / 2 + (j - 3) * dx
                if x > 0.0 #&& x <0.1
                    SigmaS[i, j] = v
                end

                #				if x>0.4 && x<0.5
                #					SigmaT[i,j] = v
                #				end
                #
                #
                #				if x>0.7 && x<0.8
                #					SigmaT[i,j] = v
                #				end
            end
        end
        SigmaT = SigmaS + SigmaA


    elseif ThisICType == "checkerboard"
        x0 = 0.0
        x1 = 7.0
        y0 = 0.0
        y1 = 7.0
        dx = (x1 - x0) / (nx)
        dy = (y1 - y0) / (ny)
        SigmaS = zeros(ny + 4, nx + 4)
        SigmaT = zeros(ny + 4, nx + 4)
        phi0 = zeros(Float64, nquadpoints, ny + 4, nx + 4)
        tEnd = 3.2
        Source = zeros(Float64, nquadpoints, ny + 4, nx + 4)

        for j = 3:nx+2
            for i = 3:ny+2
                y = y0 + dy / 2 + (i - 3) * dy
                x = x0 + dx / 2 + (j - 3) * dx
                cx = ceil(x)
                cy = ceil(y)
                g =
                    float(ceil((cx + cy) / 2.0) * 2.0 == (cx + cy)) * (
                        (1.0 < cx) &&
                        (cx < 7.0) &&
                        (1.0 < cy) &&
                        ((cy - 2.0 * abs(cx - 4.0)) < 4.0)
                    )
                SigmaS[i, j] = (1 - g) * 1 + g * 0
                SigmaT[i, j] = (1 - g) * 0 + g * 10 + SigmaS[i, j]
                if ((3 < x) && (x < 4) && (3 < y) && (y < 4))
                    for q = 1:nquadpoints
                        Source[q, i, j] = 1.0
                    end
                    SigmaS[i, j] = 0.0
                    SigmaT[i, j] = 10.0
                end
            end
        end
    elseif ThisICType == "dim1"
        SigmaS = 1 * ones(ny + 4, nx + 4)
        SigmaA = 0 * ones(ny + 4, nx + 4)
        SigmaT = SigmaS + SigmaA
        x0 = -1.5
        x1 = 1.5
        y0 = -1.5
        y1 = 1.5
        dx = (x1 - x0) / (nx)
        dy = (y1 - y0) / (ny)
        tEnd = 10.0
        phi0 = setInitialConditionsGaussian1d(nx, ny, nquadpoints, dx, dy, x0, y0)
        Source = zeros(Float64, nquadpoints, ny + 4, nx + 4)
    else
        #warning('Not yet implemented')
    end
    return x0, x1, y0, y1, tEnd, SigmaS, SigmaT, Source, phi0
end



function setInitialConditionsGaussian(nx, ny, nq, dx, dy, domainx0, domainy0)
    # init to zero
    phi = zeros(Float64, nq, ny + 4, nx + 4)
    x0 = 0.0
    y0 = 0.0
    s2 = 0.03^2
    floor = 1e-4
    f(x, y) = max(
        floor,
        1.0 / (4.0 * pi * s2) *
        exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / 4.0 / s2),
    )
    for j = 3:nx+2
        for i = 3:ny+2
            y = domainy0 + dy / 2 + (i - 3) * dy
            x = domainx0 + dx / 2 + (j - 3) * dx
            for q = 1:nq
                phi[q, i, j] = f(x, y) / 4.0 / pi
            end
        end
    end
    return phi
end

function setInitialConditionsGaussian1d(nx, ny, nq, dx, dy, domainx0, domainy0)
    # init to zero
    phi = zeros(Float64, nq, ny + 4, nx + 4)
    x0 = 0.0
    y0 = 0.0
    s2 = 0.03^2
    floor = 1e-4
    f(x, y) = max(
        floor,
        1.0 / (4.0 * pi * s2) *
        exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / 4.0 / s2),
    )
    for j = 3:nx+2
        for i = 3:ny+2
            y = domainy0 + dy / 2 + (i - 3) * dy
            x = domainx0 + dx / 2 + (j - 3) * dx
            for q = 1:nq
                phi[q, i, j] = f(x, 0) / 4.0 / pi
            end
        end
    end
    return phi
end


function setInitialConditionsBeam(
    nx,
    ny,
    nq,
    dx,
    dy,
    domainx0,
    domainy0,
    quadpoints,
    beamdirection,
)
    # init to zero
    phi = 1e-14 .* ones(Float64, nq, ny + 4, nx + 4)
    # beam postion:
    x0 = -0.50
    y0 = 0.00
    # beam width, height
    height = 0.05
    width = 0.05

    for j = 3:nx+2
        for i = 3:ny+2
            y = domainy0 + dy / 2 + (i - 3) * dy
            x = domainx0 + dx / 2 + (j - 3) * dx
            if x > x0 - width && x < x0 + width && y > y0 - height && y < y0 + height # beam should be here!
                # beam is a gaussian with pole pointing into beamdirection
                for q = 1:nq
                    Omega = quadpoints[q, :]
                    eps = 1^2
                    x = dot(Omega, beamdirection) - 1
                    val = 1.0 / sqrt(2.0 * pi * eps) * exp(-x^2 / (2 * eps))
                    phi[q, i, j] = val
                end
            end
        end
    end
    return phi
end
