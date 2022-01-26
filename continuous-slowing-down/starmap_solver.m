function [output,par] = starmap_solver(par)
%STARMAP_SOLVER
%   A second order staggered grid finite difference solver for
%   linear hyperbolic moment approximations (such as P_N and
%   SP_N) to the equations of radiative transfer in a 2D
%   geometry. Needs to be called from an example file by
%   STARMAP_SOLVER(PAR), where PAR defines the problem
%   parameters and system matrices.
%   OUTPUT = STARMAP_SOLVER(PAR) writes the final state of the
%   zeroth moment and its (x,y) coordinates into the struct
%   OUTPUT.
%
%   Version 1.5
%   Copyright (c) 05/31/2014 Benjamin Seibold and Martin Frank
%   http://www.math.temple.edu/~seibold
%   http://www.mathcces.rwth-aachen.de/5people/frank/start
%   Contributers: Edgar Olbrant (v1.0), Kerstin Kuepper (v1.5)
%
%   Date: March 2016.
%
%   StaRMAP project website:
%   http://math.temple.edu/~seibold/research/starmap/

%========================================================================
% License
%========================================================================
%   Permission is hereby granted, free of charge, to any person obtaining
%   a copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to
%   permit persons to whom the Software is furnished to do so, subject to
%   the following conditions:
%   * The above copyright notice and this permission notice shall be
%     included in all copies or substantial portions of the Software.
%   * Appropriate references shall be given to the name of the Software
%     (StaRMAP), the website of the Software
%     (http://www.math.temple.edu/~seibold/research/starmap)
%     and the research article (http://arxiv.org/abs/1211.2205)
%     corresponding to the Software.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
%   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
%   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
%   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.

%========================================================================
% Parameter Defaults
%========================================================================
tic
if ~isfield(par,'Mx')||~isfield(par,'My') ||~isfield(par,'Mz')    % If no
    par.Mx = [0 1 0 0;1 0 0 0;0 0 0 0;0 0 0 0]/sqrt(3); % system matrices
    par.n_mom = 1;                                   % specified, use P1.
    par.My = par.Mx([1 3 2 4],[1 3 2 4]);
    par.Mz = par.Mx([1 3 4 2],[1 3 4 2]);
end
if ~isfield(par,'n'),        par.n = [100 100 100]; end
% If 1D or 2D test case, modify closure matrices.
if par.n(1)==1                % Find redundancies, due to x-independence.      
    [index1,index2] = find(par.Mx);
    i1 = index1(index1>index2); i2 = index2(index1>index2);
    k = 1;
    while k<length(i1)
        temp = i1(k); i1 = i1(i2~=temp); i2 = i2(i2~=temp); k = k+1;
    end
    index = setdiff(1:(par.n_mom+1)^2,i1);
    if par.n(2)==1            % Find redundancies, due to y-independence.
        [index1,index2] = find(par.My);
        i1 = index1(index1>index2); i2 = index2(index1>index2);
        k = 1;
        while k<length(i1)
            temp = i1(k); i1 = i1(i2~=temp); i2 = i2(i2~=temp); k = k+1;
        end
        index = setdiff(index,i1);
        k = 2;               % Find redundancies, due to xy-independence.
        while k<=length(index)
            if par.Mz(index(k-1),index(k))==0, 
                index = setdiff(index,index(k));
            else k = k+1;
            end
        end
    end 
    par.Mx = par.Mx(index,index);             % Delete vanishing moments.
    par.My = par.My(index,index);
    par.Mz = par.Mz(index,index);
    if isfield(par,'mom_order')
        par.mom_order = par.mom_order(index);       % Adapt moment order.
    end
end
n_sys = length(par.Mx);                    % Number of system components.
if ~isfield(par,'name'),     par.name = 'Output'; end
if ~isfield(par,'ic'),       par.ic = @zero; end
if ~isfield(par,'sigma_a'),  par.sigma_a = @zero; end
if ~isfield(par,'sigma_s0'), par.sigma_s0 = @zero; end
flag_filter = isfield(par,'filterfunction')&&isfield(par,'f_position')...
    &&isfield(par,'mom_order');
if flag_filter, f_mom_order = par.mom_order;
else par.filterfunction = @ones; par.f_position = '';
end
flag_ani_scatter = isfield(par,'sigma_sm')&&isfield(par,'mom_order');
if ~flag_ani_scatter           % If anisotropic part of scattering kernel
    par.mom_order = [0 ones(1,n_sys-1)]; par.sigma_sm = @zero;      % not
end                    % available, set up data for isotropic scattering.
if ~isfield(par,'source'),  par.source = @zero; end
if ~isfield(par,'source_ind') % If no source moment vector specified, use
    par.source_ind = 1:((nargin(par.source)>=0)+...   % zeroth moment for
        (nargin(par.source)>4)*(n_sys-1));   % isotropic, and all moments
end                                             % for anisotropic source.
if ~isfield(par,'ax'),       par.ax = [0 1 0 1 0 1]; end
if ~isfield(par,'cn'),       par.cn = 0.99; end     % Default CFL number.
if ~isfield(par,'bc'),       par.bc = [0 0]; end
if ~isfield(par,'t_plot'),   par.t_plot = 0:.1:1; end
if ~isfield(par,'tfinal'),   par.tfinal = par.t_plot(end); end
if ~isfield(par,'output'),   par.output = @default_output; end
if ~isfield(par,'mom_output'), par.mom_output = 1; end % Moments plotted.
if strcmp(par.mom_output,'all'), par.mom_output = 1:n_sys; end
if ~isfield(par,'int_weight'), par.int_weight = @(m,t)0; end

%========================================================================
% Initialization
%========================================================================
% Create index matrices and vectors for grid type
Ix = cellfun(@find,num2cell(par.Mx',1),'Un',0);         % Grid dependency
Iy = cellfun(@find,num2cell(par.My',1),'Un',0);                % indices.
Iz = cellfun(@find,num2cell(par.Mz',1),'Un',0);
c111 = 1; c122 = []; c121 = []; c112 = [];        % First component is on 
c211 = []; c222 = []; c221 = []; c212 = [];                  %(1,1)-grid.
while length([c111,c122,c121,c112,c211,c222,c221,c212])<n_sys % Determine
    c211 = unique(vertcat(c211',Ix{c111},Iy{c221},Iz{c212})');  % indices
    c121 = unique(vertcat(c121',Ix{c221},Iy{c111},Iz{c122})');  % of grid
    c112 = unique(vertcat(c112',Ix{c212},Iy{c122},Iz{c111})');    % types
    c221 = unique(vertcat(c221',Ix{c121},Iy{c211},Iz{c222})');     % from
    c212 = unique(vertcat(c212',Ix{c112},Iy{c222},Iz{c211})');   % system
    c122 = unique(vertcat(c122',Ix{c222},Iy{c112},Iz{c121})');% matrices.
    c111 = unique(vertcat(c111',Ix{c211},Iy{c121},Iz{c112})');
    c222 = unique(vertcat(c222',Ix{c122},Iy{c212},Iz{c221})');
end
if isempty(c221), c221 = []; end
if isempty(c222), c222 = []; end

gtx = ones(1,n_sys); gty = gtx; gtz = gtx;                  % Grid types
gtx([c211,c212,c221,c222]) = 2;                              % in x,y,z.
gty([c121,c122,c221,c222]) = 2;
gtz([c112,c122,c212,c222]) = 2;

% Grid size and time step
h = (par.ax([2 4 6])-par.ax([1 3 5]))./par.n;                % Grid size.
dt = min(h)/abs(eigs(par.Mz,1,'lm'))*par.cn;         % Maximal time step.
if par.n(1)~=1 % Adapt dt to 3D test case.
    dt = dt/sqrt(3);
elseif par.n(2)~=1 % Adapt dt to 2D test case.
    dt = dt/sqrt(2);
end


% Create staggered grids and coefficients
x{1} = par.ax(1)+h(1)/2:h(1):par.ax(2)-h(1)/2;              % Cell center
y{1} = par.ax(3)+h(2)/2:h(2):par.ax(4)-h(2)/2;             % coordinates.
z{1} = par.ax(5)+h(3)/2:h(3):par.ax(6)-h(3)/2;
x{2} = par.ax(1)+h(1)*(1-par.bc(1)):h(1):par.ax(2);           % Staggered
y{2} = par.ax(3)+h(2)*(1-par.bc(2)):h(2):par.ax(4);        % coordinates.
z{2} = par.ax(5)+h(3)*(1-par.bc(3)):h(3):par.ax(6);

X = cell(2,2,2); Y = X; Z = X;                         % Staggered grids.
[X(:,:,1),Y(:,:,1),Z(:,:,1)] = cellfun(@ndgrid,[x;x]',[y;y],...
    {z{1},z{1};z{1},z{1}},'Un',0);
[X(:,:,2),Y(:,:,2),Z(:,:,2)] = cellfun(@ndgrid,[x;x]',[y;y],...
    {z{2},z{2};z{2},z{2}},'Un',0);

n1 = par.n; n2 = size(X{2,2,2});              % Sizes of staggered grids.
extendx = {[1:par.bc(1),1:n1(1),par.bc(1)*(n1(1)-1)+1],...    % Extension
    [n2(1)*(1:1-par.bc(1)),1:n2(1)]}; % indices for
extendy = {[1:par.bc(2),1:n1(2),par.bc(2)*(n1(2)-1)+1],...     % boundary
    [n2(2)*(1:1-par.bc(2)),1:n2(2)]};                       % conditions.
extendz = {[1:par.bc(3),1:n1(3),par.bc(3)*(n1(3)-1)+1],...
    [n2(3)*(1:1-par.bc(3)),1:n2(3)]};
sg = [gtx;gty;gtz;par.mom_order]';             % Construct staggered grid
sg = unique(sg(2:end,:),'rows')';            % indices with moment order.

% Initialize field variables
sS = cell(2,2,2,max(par.mom_order));                 % Scattering fields.
sT = sS; ET = sS;
U = cell(1,n_sys); dxU = U; dyU = U; dzU = U; Q = U;          % Unknowns.
for j = 1:n_sys                           % If initial condition only two 
    U{j} = X{gtx(j),gty(j),gtz(j)}*0+...                 % arguments, set
        capargs(par.ic,X{gtx(j),gty(j),gtz(j)},...   % all higher moments 
        Y{gtx(j),gty(j),gtz(j)},Z{gtx(j),gty(j),gtz(j)},j)*... % to zero.
        (nargin(par.ic)>3||j==1);           
    Q{j} = X{gtx(j),gty(j),gtz(j)}*0;%Initialize source terms with zeros.
end

% Initialize filtering variables
if flag_filter
    if strcmp(par.f_position,'full')
        filterTerm = cell(2,2,2,max(f_mom_order));        % Filter term.
    elseif strcmp(par.f_position,'substep')
        filterTerm = cell(2,2,2,max(f_mom_order),2);
    end
    f_sg = [gtx;gty;gtz;f_mom_order]';     % Staggered grid indices with 
    f_sg = unique(f_sg(2:end,:),'rows')';                % moment order.
end
% Initialize integral with zeros.
Int = Q;

% Initialize density function.
if ~isfield(par,'density') % If no density specified set to 1.
    Rho{1,1,1} = 1; Rho{1,2,1} = 1; Rho{1,1,2} = 1; Rho{1,2,2} = 1;
    Rho{2,1,1} = 1; Rho{2,2,1} = 1; Rho{2,1,2} = 1; Rho{2,2,2} = 1;
    density_extendx = {1 1}; density_extendy = {1 1}; density_extendz = {1 1};
else
    Rho = cell(2,2,2); rho_min = 1;
    for i = 1:2
        for j = 1:2
            for k = 1:2
                Rho{i,j,k} = X{i,j,k}*0 + par.density(X{i,j,k},Y{i,j,k},Z{i,j,k},par);
                rho_min = min(rho_min,min(min(min((abs(Rho{i,j,k})))))); % Minimal density.
            end
        end
    end
    density_extendx = extendx; density_extendy = extendy; density_extendz = extendz;
    % Adjust time step        
    dt = rho_min*dt;
end

dt = par.tfinal/ceil(par.tfinal/dt);    % Ensure integer number of steps.

fprintf(['CPU-time Initialization: ',num2str(toc),'s. \n'])
%========================================================================
% Time Loop
%========================================================================
par.t_plot = [par.t_plot inf]; plot_count = 1;
cputime = zeros(1,4);
for t = dt:dt:par.tfinal
    U0 = U(par.mom_output);                     % Store current solution.
    for j = par.source_ind     % Source: evaluate only active components.
        Q{j} = capargs(par.source,X{gtx(j),gty(j),gtz(j)},...  % Evaluate 
            Y{gtx(j),gty(j),gtz(j)},Z{gtx(j),gty(j),gtz(j)},...  % source
            t-dt/2,j);                            % at half-time of step.
    end
    % Evaluate material parameters
    tic
    evaluate = (t==dt|[nargin(par.sigma_a)>3,nargin(par.sigma_s0)>3,...
        nargin(par.sigma_sm)>4])&[1 1 flag_ani_scatter];
    if evaluate(1)                                % Absorption: evaluate, 
        sA = cellfun(@(x,y,z)capargs(par.sigma_a,...    % if time-dependent
            x,y,z,t-dt/2),X,Y,Z,'Un',0);            % or first time step.
    end
    if evaluate(2)                     % Scattering, isotropic component:
        for j = 1:size(sg,2)       % evaluate, if time-dependent or first
            sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = capargs(par.sigma_s0,... % time
                X{sg(1,j),sg(2,j),sg(3,j)},Y{sg(1,j),sg(2,j),sg(3,j)},...
                Z{sg(1,j),sg(2,j),sg(3,j)},t-dt/2);               % step.
        end
    end
    if evaluate(3)                  % Scattering, anisotropic components:
        for j = 1:size(sg,2)             % evaluate, if time-dependent or
            sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...  % first time step.
                sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)}-capargs(par.sigma_sm,...
                X{sg(1,j),sg(2,j),sg(3,j)},Y{sg(1,j),sg(2,j),sg(3,j)},...
                Z{sg(1,j),sg(2,j),sg(3,j)},sg(4,j),t-dt/2);
        end
    end
    if any(evaluate)     % Total sigma_t^m = sigma_a+sigma_s^0-sigma_s^m.
        for j = 1:size(sg,2)                 % (Re)compute, if any of the 
            sT{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...             % above
                sA{sg(1,j),sg(2,j),sg(3,j)}+...           %(re)evaluated.
                sS{sg(1,j),sg(2,j),sg(3,j),sg(4,j)};
            ET{sg(1,j),sg(2,j),sg(3,j),sg(4,j)} = ...% Total sigma decay.
                expm1div(-sT{sg(1,j),sg(2,j),sg(3,j),sg(4,j)}*dt/2);
        end
        EA = expm1div(-sA{1,1}*dt/2);     % Decay term for zeroth moment.
    end
    % Evaluate filter function
    if (t==dt||nargin(par.filterfunction)>6)&&flag_filter  % Evaluate, if
        for j = 1:size(f_sg,2)       % first time step or time-dependent.
            if strcmp(par.f_position,'full')
                filterTerm{f_sg(1,j),f_sg(2,j),f_sg(3,j),f_sg(4,j)} = ...
                    capargs(par.filterfunction,par,dt,f_sg(4,j),...
                    X{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Y{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Z{f_sg(1,j),f_sg(2,j),f_sg(3,j)},t);
            elseif strcmp(par.f_position,'substep')
                filterTerm{f_sg(1,j),f_sg(2,j),f_sg(3,j),f_sg(4,j),1} = ...
                    capargs(par.filterfunction,par,dt,f_sg(4,j),...
                    X{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Y{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Z{f_sg(1,j),f_sg(2,j),f_sg(3,j)},t-dt/2);
                filterTerm{f_sg(1,j),f_sg(2,j),f_sg(3,j),f_sg(4,j),2} = ...
                    capargs(par.filterfunction,par,dt,f_sg(4,j),...
                    X{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Y{f_sg(1,j),f_sg(2,j),f_sg(3,j)},...
                    Z{f_sg(1,j),f_sg(2,j),f_sg(3,j)},t);
            end
        end
    end
    cputime(2) = cputime(2)+toc;
    % Update solution
    k=1;
    for step = [1 2 1]                   % Strang splitting of sub-steps.
        switch step
        case 1, tic    % Half step with (111), (221), (212), (122) grids.
        for j = c111                         % Compute update from (111).
            dxU{j} = diff(U{j}(extendx{1},:,:)./...
                Rho{1,1,1}(density_extendx{1},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{1},:)./...
                Rho{1,1,1}(:,density_extendy{1},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{1})./...
                Rho{1,1,1}(:,:,density_extendz{1}),[],3)/h(3);
        end
        for j = c221                         % Compute update from (221).
             dxU{j} = diff(U{j}(extendx{2},:,:)./...
                Rho{2,2,1}(density_extendx{2},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{2},:)./...
                Rho{2,2,1}(:,density_extendy{2},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{1})./...
                Rho{2,2,1}(:,:,density_extendz{1}),[],3)/h(3);
        end
        for j = c212                         % Compute update from (212).
            dxU{j} = diff(U{j}(extendx{2},:,:)./...
                Rho{2,1,2}(density_extendx{2},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{1},:)./...
                Rho{2,1,2}(:,density_extendy{1},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{2})./...
                Rho{2,1,2}(:,:,density_extendz{2}),[],3)/h(3);
        end
        for j = c122                         % Compute update from (122).
            dxU{j} = diff(U{j}(extendx{1},:,:)./...
                Rho{1,2,2}(density_extendx{1},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{2},:)./...
                Rho{1,2,2}(:,density_extendy{2},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{2})./...
                Rho{1,2,2}(:,:,density_extendz{2}),[],3)/h(3);
        end
        for j = [c211 c121 c112 c222]       % Update components on grids.
            W = -sumcell([dxU(Ix{j}),dyU(Iy{j}),dzU(Iz{j})],...
                [par.Mx(j,Ix{j}),par.My(j,Iy{j}),par.Mz(j,Iz{j})]);
            U{j} = U{j}+dt/2*(W+Q{j}-...
                sT{gtx(j),gty(j),gtz(j),par.mom_order(j)}.*U{j}).*...
                ET{gtx(j),gty(j),gtz(j),par.mom_order(j)};
            if strcmp(par.f_position,'substep')       % Multiplication by
                U{j} = U{j}.*filterTerm{gtx(j),gty(j),gtz(j),f_mom_order(j),k};
            end                                    % substep filter term.
        end
        cputime(1) = cputime(1)+toc;
        case 2, tic    % Half step with (211), (121), (112), (222) grids.
        for j = c222                         % Compute update from (222).
            dxU{j} = diff(U{j}(extendx{2},:,:)./...
                Rho{2,2,2}(density_extendx{2},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{2},:)./...
                Rho{2,2,2}(:,density_extendy{2},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{2})./...
                Rho{2,2,2}(:,:,density_extendz{2}),[],3)/h(3);
        end
        for j = c112                         % Compute update from (112).
            dxU{j} = diff(U{j}(extendx{1},:,:)./...
                Rho{1,1,2}(density_extendx{1},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{1},:)./...
                Rho{1,1,2}(:,density_extendy{1},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{2})./...
                Rho{1,1,2}(:,:,density_extendz{2}),[],3)/h(3);
        end
        for j = c121                         % Compute update from (121).
            dxU{j} = diff(U{j}(extendx{1},:,:)./...
                Rho{1,2,1}(density_extendx{1},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{2},:)./...
                Rho{1,2,1}(:,density_extendy{2},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{1})./...
                Rho{1,2,1}(:,:,density_extendz{1}),[],3)/h(3);
        end
        for j = c211                         % Compute update from (211).
            dxU{j} = diff(U{j}(extendx{2},:,:)./...
                Rho{2,1,1}(density_extendx{2},:,:),[],1)/h(1);
            dyU{j} = diff(U{j}(:,extendy{1},:)./...
                Rho{2,1,1}(:,density_extendy{1},:),[],2)/h(2);
            dzU{j} = diff(U{j}(:,:,extendz{1})./...
                Rho{2,1,1}(:,:,density_extendz{1}),[],3)/h(3);
        end
        for j = [c111 c221 c212 c122]       % Update components on grids.
            W = -sumcell([dxU(Ix{j}),dyU(Iy{j}),dzU(Iz{j})],...
                [par.Mx(j,Ix{j}),par.My(j,Iy{j}),par.Mz(j,Iz{j})]);
            for k = 1:2                         % Perform two half-steps.
                if j==1        % Zeroth moment decays by absorption only.
                    U{j} = U{j}+dt/2*(W+Q{j}-sA{1,1,1}.*U{j}).*EA;
                else                  % All other moments decay normally.
                    U{j} = U{j}+dt/2*(W+Q{j}-...
                        sT{gtx(j),gty(j),gtz(j),par.mom_order(j)}.*...
                        U{j}).*ET{gtx(j),gty(j),gtz(j),par.mom_order(j)};
                    if strcmp(par.f_position,'substep')  % Multiplication
                        U{j} = U{j}.*...        % by substep filter term.
                            filterTerm{gtx(j),gty(j),gtz(j),f_mom_order(j),k};
                    end
                end
            end
        end
        cputime(1) = cputime(1)+toc;
        end
    end
    if strcmp(par.f_position,'full') % Multiplication by full filter term 
        for j= 2:n_sys                         % at the end of time step.
            U{j} = U{j}.*filterTerm{gtx(j),gty(j),gtz(j),f_mom_order(j)};
        end
    end
    % Compute weighted time integral using the trapezoidal rule.
    tic
    weight = par.int_weight(1:n_sys,t-dt/2)*dt/(1+(t==dt||t==par.tfinal));
    for j = 1:n_sys
        Int{j} = Int{j}+U{j}.*weight(j);
    end
    cputime(4) = cputime(4)+toc;
    % Plotting
    tic
    while t>=par.t_plot(plot_count)-1e-14  % If current time has exceeded
        lambda = (par.t_plot(plot_count)-t+dt)/dt; %plotting time, define
        Uplot = cellfun(@(x,y)(1-lambda)*x+lambda*y,...    % solution via
            U0,U(par.mom_output),'Un',0); % linear interpolation in time.
        xplot = x(gtx(par.mom_output));                 % Assign grids at
        yplot = y(gty(par.mom_output));              % outputted moments.
        zplot = z(gty(par.mom_output));
        if length(par.mom_output)==1   % If only a single moment plotted,
            xplot = xplot{:}; yplot = yplot{:}; zplot = zplot{:};...
            Uplot = Uplot{:};                    % remove cell structure.
        end                                             
        if nargout(par.output)                     % Call output routine,
            par = par.output(par,xplot,yplot,zplot,Uplot,plot_count);
        else                       % allowing for it to modify the struct
            par.output(par,xplot,yplot,zplot,Uplot,plot_count)     % par.
        end
        plot_count = plot_count+1;
        fprintf('CPU-time:%12.2fs\n',sum(cputime([1,2,4])))
        fprintf('Remaining:%11.2fs\n',sum(cputime([1,2,4])*(par.tfinal./t-1)))
    end
    cputime(3) = cputime(3)+toc;
    % fprintf(['t = ',num2str(t),' CPU-time %12.2fs %12.2fs %12.2fs %12.2fs \n'],cputime)
end
fprintf('%s with %s%d\n %d moments, %dx%dx%d grid, %0.0f time steps\n',... 
    par.name,par.closure,par.n_mom,n_sys,par.n,par.tfinal/dt)% Print test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%12.2fs%5.0f%%\n ',...   % and CPU times.
    'absorb/scatter:%7.2fs%5.0f%%\n ','plotting:%13.2fs%5.0f%%\n ',...
    'integration:%10.2fs%5.0f%%\n'],cputime)

% Output final state of solution
if nargout, output = struct('x',x(gtx),'y',y(gty),'z',z(gtz),'U',U,'Int',Int); end

%========================================================================
% Technical Functions
%========================================================================
function S = sumcell(A,w)
% Add vector of cells A, weighted with vector w.
S = A{1}*w(1); for j = 2:length(w), S = S+A{j}*w(j); end

function z = capargs(fct,varargin)
% Call function fct with as many arguments as it requires (at least 1),
% and ignore further arguments.
narg = max(nargin(fct),1);
z = fct(varargin{1:narg});

function f = zero(varargin)
% Zero function.
f = zeros(size(varargin{1}));

function y = expm1div(x)
% Function (exp(x)-1)/x that is accurate for x close to zero.
y = 1+x*.5+x.^2/6;
ind = abs(x)>2e-4;
y(ind) = (exp(x(ind))-1)./x(ind);

function default_output(par,x,y,z,U,step)
% Default plotting routine.
clf, imagesc(x,y,U(:,:,ceil(length(z)/2))'), axis xy equal tight
title(sprintf('%s at t = %0.2f',par.name,par.t_plot(step))), drawnow
