# this file is not used in the code. it contains all the testcases but we now use getProblemSpecificParameters instead
# you can however move a testcase from here to the getProblemSpecificParameters file but have to remove the KS. everywhere




function setTestCaseInitialCondition!(KS::KineticSolver,phi::Array{Float64,3})
    
    if KS.ictype=="linesource"
        KS.SigmaS = 1*ones(KS.ny,KS.nx);
        SigmaA = 0*ones(KS.ny,KS.nx);
        KS.SigmaT = KS.SigmaS+SigmaA;
        KS.x0 = -1.5
        KS.x1 =  1.5
        KS.y0 = -1.5
        KS.y1 =  1.5
      
        setInitialConditionsGaussian(KS,phi)
    elseif KS.ictype=="checkerboard"
        KS.x0 = 0.0;
        KS.x1 = 7.0;
        KS.y0 = 0.0;
        KS.y1 = 7.0;
        KS.dx = (KS.x1-KS.x0)/(KS.nx-4)
        KS.dy = (KS.y1-KS.y0)/(KS.ny-4)
        KS.SigmaS = zeros(KS.ny,KS.nx);
        KS.SigmaT = zeros(KS.ny,KS.nx);
        phi[:] = 1e-10;
        KS.tEnd = 3.2
        KS.dt = KS.cfl/2*(KS.dx*KS.dy)/(KS.dx+KS.dy)

        KS.nt = ceil(KS.tEnd/KS.dt);


       # KS.rangex = collect(1:KS.nx)+2 # two ghost cells left and right
       # KS.rangey = collect(1:KS.ny)+2 # same
        for j=KS.rangex
            for i=KS.rangey   
                y = KS.y0+KS.dy/2 + (i-3)*KS.dy
                x = KS.x0+KS.dx/2 + (j-3)*KS.dx     
                cx = ceil(x);
                cy = ceil(y);
                #cx = x;
                #cy = y;
                g = float(ceil((cx+cy)/2.)*2.==(cx+cy))*((1.<cx)&&(cx<7.)&&(1.<cy)&&((cy-2.*abs(cx-4.0))<4.0));
                KS.SigmaS[i,j] = (1-g)*1+g*0;   
                KS.SigmaT[i,j] = (1-g)*0+g*10 + KS.SigmaS[i,j]
                
                if ((3<x)&&(x<4)&&(3<y)&&(y<4))
                        for q=KS.rangequad
                                  KS.Source[q,i,j] = 1.0;#*KS.dx*KS.dy;#*4.0*pi/KS.nquadpoints;
                          end
                    end 
            end
        end
    elseif KS.ictype=="checkerboard2"
        KS.x0 = 0.0;
        KS.x1 = 7.0;
        KS.y0 = 0.0;
        KS.y1 = 7.0;
        KS.dx = (KS.x1-KS.x0)/(KS.nx-4)
        KS.dy = (KS.y1-KS.y0)/(KS.ny-4)
        KS.SigmaS = zeros(KS.ny,KS.nx);
        KS.SigmaT = zeros(KS.ny,KS.nx);
        phi[:] = 1e-12;
        KS.tEnd = 3.2
        KS.dt = KS.cfl/2*(KS.dx*KS.dy)/(KS.dx+KS.dy)

        KS.nt = ceil(KS.tEnd/KS.dt);


       # KS.rangex = collect(1:KS.nx)+2 # two ghost cells left and right
       # KS.rangey = collect(1:KS.ny)+2 # same
        for j=KS.rangex
            for i=KS.rangey   
                y = KS.y0+KS.dy/2 + (i-3)*KS.dy
                x = KS.x0+KS.dx/2 + (j-3)*KS.dx     
                cx = ceil(x);
                cy = ceil(y);
                #cx = x;
                #cy = y;
                g = float(ceil((cx+cy)/2.)*2.==(cx+cy))*((1.<cx)&&(cx<7.)&&(1.<cy)&&((cy-2.*abs(cx-4.0))<4.0));
                KS.SigmaS[i,j] = (1-g)*1+g*0;   
                KS.SigmaT[i,j] = (1-g)*0+g*10 + KS.SigmaS[i,j]
                
                if ((3<x)&&(x<4)&&(3<y)&&(y<4))
                        for q=KS.rangequad
                                  KS.Source[q,i,j] = 1.0;#KS.nquadpoints;#*KS.dx*KS.dy;#*4.0*pi/KS.nquadpoints;
                        end
                        KS.SigmaS[i,j] = 0.0;
                        KS.SigmaT[i,j] = 10.0;
                    end 
            end
        end
    elseif KS.ictype=="beam"
        KS.x0 = 0.0;
        KS.x1 = 7.0;
        KS.y0 = 0.0;
        KS.y1 = 7.0;
        KS.dx = (KS.x1-KS.x0)/(KS.nx-4)
        KS.dy = (KS.y1-KS.y0)/(KS.ny-4)
        KS.SigmaS = zeros(KS.ny,KS.nx);
        KS.SigmaT = zeros(KS.ny,KS.nx);
        phi[:] = 1e-8;
        KS.tEnd = 4.0
        KS.dt = KS.cfl/2*(KS.dx*KS.dy)/(KS.dx+KS.dy)

        KS.nt = ceil(KS.tEnd/KS.dt);
        
        tmpMax = -10.0
        tmpMid = -5.0;
        qMax = -1;
        qMid = -1;
        for q=KS.rangequad
            if KS.Q.pointsxyz[q,1]>= tmpMax
                tmpMid = tmpMax;
                tmpMax = KS.Q.pointsxyz[q,1]; 
                qMid = qMax;
                qMax = q;
            elseif KS.Q.pointsxyz[q,1]>= tmpMid
                tmpMid = KS.Q.pointsxyz[q,1]; 
                qMid = q;
            end
        end
        #println(KS.Q.pointsxyz[qMid,1]) 
        #println(KS.Q.pointsxyz[qMax,1]) 
        #error("Max directions are ",qMax, " and ", qMid);


       # KS.rangex = collect(1:KS.nx)+2 # two ghost cells left and right
       # KS.rangey = collect(1:KS.ny)+2 # same
        for j=KS.rangex
            for i=KS.rangey   
                y = KS.y0+KS.dy/2 + (i-3)*KS.dy
                x = KS.x0+KS.dx/2 + (j-3)*KS.dx     
                cx = ceil(x);
                cy = ceil(y);
                #cx = x;
                #cy = y;
                #g = float(ceil((cx+cy)/2.)*2.==(cx+cy))*((1.<cx)&&(cx<7.)&&(1.<cy)&&((cy-2.*abs(cx-4.0))<4.0));
                KS.SigmaS[i,j] = 0.0#(1-g)*1+g*0;   
                KS.SigmaT[i,j] = 0.1#(1-g)*0+g*10 + KS.SigmaS[i,j]
                

                if (3.45<x && x< 3.55 &&i==20)
                    for q=KS.rangequad
                        #if q==qMid || q == qMax
                        if KS.Q.pointsxyz[q,1] > 0.7
                        #if KS.Q.pointsxyz[q,1]>=0
                        #KS.Source[q,i,j] = 10.0;#KS.nquadpoints;#*KS.dx*KS.dy;#*4.0*pi/KS.nquadpoints;
                        phi[q,i,j] = 10.0;
                        end
                    end
                end 
            end
        end
    elseif KS.ictype=="circle"
        KS.x0 = -4.0;
        KS.x1 = 4.0;
        KS.y0 = -4.0;
        KS.y1 = 4.0;
        KS.dx = (KS.x1-KS.x0)/(KS.nx-4)
        KS.dy = (KS.y1-KS.y0)/(KS.ny-4)
        KS.SigmaS = zeros(KS.ny,KS.nx);
        KS.SigmaT = zeros(KS.ny,KS.nx);
        phi[:] = 1e-5;
        KS.tEnd = 40
        KS.dt = KS.cfl/2*(KS.dx*KS.dy)/(KS.dx+KS.dy)

        KS.nt = ceil(KS.tEnd/KS.dt);


       # KS.rangex = collect(1:KS.nx)+2 # two ghost cells left and right
       # KS.rangey = collect(1:KS.ny)+2 # same
        for j=KS.rangex
            for i=KS.rangey   
                y = KS.y0+KS.dy/2 + (i-3)*KS.dy
                x = KS.x0+KS.dx/2 + (j-3)*KS.dx     
                cx = ceil(x);
                cy = ceil(y);
                #cx = x;
                #cy = y;
                KS.SigmaS[i,j] = 1# 1.81-0.09056
                KS.SigmaT[i,j] = 1.5#1.81
                r = 1.0;
                dr = 0.2; # width of the circle
                if (x*x+y*y)>(r-dr) && (x*x+y*y)<(r+dr)
                        for q=KS.rangequad
                                  KS.Source[q,i,j] = 1.0;#KS.nquadpoints;#*KS.dx*KS.dy;#*4.0*pi/KS.nquadpoints;
                        end
                    end 
            end
        end

    elseif KS.ictype=="spiral"
        KS.x0 = -4.0;
        KS.x1 = 4.0;
        KS.y0 = -4.0;
        KS.y1 = 4.0;
        KS.dx = (KS.x1-KS.x0)/(KS.nx-4)
        KS.dy = (KS.y1-KS.y0)/(KS.ny-4)
        KS.SigmaS = zeros(KS.ny,KS.nx);
        KS.SigmaT = zeros(KS.ny,KS.nx);
        phi[:] = 1e-5;
        KS.tEnd = 4
        KS.dt = KS.cfl/2*(KS.dx*KS.dy)/(KS.dx+KS.dy)

        KS.nt = ceil(KS.tEnd/KS.dt);


       # KS.rangex = collect(1:KS.nx)+2 # two ghost cells left and right
       # KS.rangey = collect(1:KS.ny)+2 # same
       
        for j=KS.rangex
            for i=KS.rangey   
                y = KS.y0+KS.dy/2 + (i-3)*KS.dy
                x = KS.x0+KS.dx/2 + (j-3)*KS.dx     
                cx = ceil(x);
                cy = ceil(y);
                #cx = x;
                #cy = y;.
                KS.SigmaS[i,j] = 1  
                KS.SigmaT[i,j] = 0
               
               # if sqrt(x*x+y*y)>sqrt(theta)-dphi && sqrt(x*x+y*y)<sqrt(theta)+dphi;
               #     for q=KS.rangequad
               #       KS.Source[q,i,j] = 1.0;#KS.nquadpoints;#*KS.dx*KS.dy;#*4.0*pi/KS.nquadpoints;
               #     end
               # end 

            end
        end

        for theta = 0:0.01:3*pi
            for r = theta/3:0.01:theta/3+0.1;
            for j=KS.rangex
                for i=KS.rangey   
                    y = KS.y0+KS.dy/2 + (i-3)*KS.dy
                    x = KS.x0+KS.dx/2 + (j-3)*KS.dx     

                    xx = cos(theta)*r/2
                    yy = sin(theta)*r/2

                    inside = (x<=xx && xx<=x+KS.dx) && (y<=yy && yy<=y+KS.dy) 
                    if inside
                        for q=KS.rangequad
                          KS.Source[q,i,j] = 1.0;#KS.nquadpoints;#*KS.dx*KS.dy;#*4.0*pi/KS.nquadpoints;
                        end
                    end 
                end
            end
            end
        end

    end


    nicewrite("SigmaS.txt",KS.SigmaS)
    nicewrite("SigmaT.txt",KS.SigmaT)

end