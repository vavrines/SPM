function [Mx,My,Mz,moment_order] = starmap_closure_spn(n_mom)
%STARMAP_CLOSURE_SPN
%   Creates SP_N moment system matrices, to be used by
%   STARMAP_SOLVER, a second order staggered grid finite
%   difference solver for linear hyperbolic moment
%   approximations to radiative transfer in 2D geometry.
%
%   The time-dependent SPN equations encoded here are
%   the ones given in
%   [Olbrant, Larsen, Frank, Seibold, JCP 238 (2013) 315-336].
%
%   Version 1.5
%   Copyright (c) 05/31/2014 Benjamin Seibold and Martin Frank
%   http://www.math.temple.edu/~seibold
%   http://www.mathcces.rwth-aachen.de/5people/frank/start
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5)
%
%   StaRMAP project website:
%   http://math.temple.edu/~seibold/research/starmap/

%   For license, see file starmap_solver.m, as published on
%   http://www.math.temple.edu/~seibold/research/starmap

%========================================================================
nstar = ceil((n_mom+1)/2);
i = 0:nstar-1;
k = 2*i.*(2*i-1)./((4*i+1).*(4*i-1));
l = 4*i.^2./((4*i+1).*(4*i-1))+...
    (2*i+1).^2./((4*i+1).*(4*i+3)).*(i<nstar-1|floor(n_mom/2)*2<n_mom);
m = 2*(2*i+1).*(i+1)./((4*i+1).*(4*i+3));
% Create matrix for x-direction
Mx = zeros(nstar*4);
Mx(nstar+(1:nstar)*3-2,1:nstar) = ...
    diag(k(2:end),-1)+diag(l)+diag(m(1:end-1),1);
Mx(1:nstar,nstar+(1:nstar)*3-2) = eye(nstar);
% Create matrix for y-direction
i = 1:nstar*3; i = [1:nstar,nstar+i+mod(i+1,3)-1];
My = Mx(i,i);
% Create matrix for z-direction
i = 1:nstar*3; i = [1:nstar,nstar+i+2*(1-mod(i+2,3))];
Mz = Mx(i,i);
% Order of moments
moment_order = [0:2:n_mom reshape([1;1;1]*((0:2:n_mom)+1),1,[])];
