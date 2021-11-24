function visualizeLattic
close all


M = NaN*ones(8,8);

val = 1.;
M(2:2:6,2) = val;
M(3:2:5,3) = val;
M(4:2:6,4) = val;
M(3:2:5,5) = val;
M(2:2:6,6) = val;
M(4,4) = 2*val;

figure()
[Ny, Nx] = size(M);
X = linspace(0, 7, Nx); Y = linspace(0, 7, Ny);

%% Figure 1: IMAGESC
figure()
h = pcolor(X, Y,M);
colormap jet
axis equal
axis tight
set(h, 'EdgeColor', 'none');
xlabel('$x$ in cm','interpreter','latex')
ylabel('$y$ in cm','interpreter','latex')
grid on
set ( gca, 'ydir', 'reverse' )
title('Layout of the lattice problem','interpreter','latex')
% yticks([7:-1:0])





end