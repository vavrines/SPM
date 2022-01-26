function [Mx,My,Mz,moment_order] = starmap_closure_pn(n_mom)
%STARMAP_CLOSURE_PN
%   Creates P_N moment system matrices, to be used by
%   STARMAP_SOLVER, a second order staggered grid finite
%   difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 2D geometry.
%
%   Version 1.1
%   Copyright (c) 07/14/2013 Benjamin Seibold and Martin Frank
%   http://www.math.temple.edu/~seibold
%   http://www.mathcces.rwth-aachen.de/5people/frank/start

%   For license, see file starmap_solver.m, as published on
%   http://www.math.temple.edu/~seibold/research/starmap

%========================================================================
n_sys = (n_mom+1)^2; % number of system components in 3D
Mx = sparse(zeros(n_sys)); My = Mx; Mz = Mx; s = size(Mx);
for m = 1:n_mom
    i = 1:2*m-1; p = (m-1)^2+i;
    v = F(m,m+1-ceil(i/2));
    Mx(sub2ind(s,p,p+2*m-1)) = v;
    My(sub2ind(s,p,p+2*m-1-(-1).^i)) = -(-1).^i.*v;
    v = B(m,m-ceil(i/2));
    Mz(sub2ind(s,p,p+2*m+1)) = v;
    i = 1:2*m-3; p = (m-1)^2+i;
    v = D(m,m-1-ceil(i/2));
    Mx(sub2ind(s,p,p+2*m+3)) = -v;
    if m>1, p(end) = p(end)+1; i(end) = i(end)+1; end
    My(sub2ind(s,p,p+2*m+3-(-1).^i)) = -(-1).^i.*v;
end
m = 1:n_mom; i = m.^2;
Mx(i,:) = sqrt(2)*Mx(i,:); My(i,:) = sqrt(2)*My(i,:);
m = 2:n_mom; i = (m+1).^2;
Mx(:,i) = sqrt(2)*Mx(:,i); My(:,i) = sqrt(2)*My(:,i);
Mx = full((Mx+Mx')/2); My = full((My+My')/2); % symmetry and factor 1/2
Mz = full(Mz+Mz');
% Order of moments in 3D
moment_order = ceil(sqrt((1:n_sys))-1);

%========================================================================

function y = B(l,m)
y = sqrt((l-m).*(l+m)/(2*l+1)/(2*l-1));

function y = D(l,m)
y = sqrt((l-m).*(l-m-1)/(2*l+1)/(2*l-1));

function y = F(l,m)
y = sqrt((l+m).*(l+m-1)/(2*l+1)/(2*l-1));
