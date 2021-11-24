function rhoNumerical = visualizeLinesource(folder, plottingparams, problemparams)
close all
FS = plottingparams.FS; % font size
LW = plottingparams.LW; % line width
PP = plottingparams.PP; % paper position

rhoExact = dlmread('exactLineSource.txt');
rhoNumerical = dlmread(strcat(folder, '/data/rhofinal.txt'))+1e-10;

rhoNumerical = rhoNumerical(3:end-2, 3:end-2); % remove ghost cells
rhoNumerical(isnan(rhoNumerical)) = 0;

[Ny, Nx] = size(rhoNumerical);
X = linspace(-1.5, 1.5, Nx); Y = linspace(-1.5, 1.5, Ny);

%% Figure 1: IMAGESC
figure()
rhoNumerical(rhoNumerical < 0) = NaN;
h = pcolor(X, Y, rhoNumerical);
set(h, 'EdgeColor', 'none');
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
xticks([]);
yticks([]);
% xticks([-1.5:0.5:1.5]); 
% yticks([-1.5:0.5:1.5]);
% grid on

hold on
line([0, 1.5; 0, 1.5/sqrt(2)]', [0, 0;0, 1.5/sqrt(2)]',  'LineWidth', 2*LW); %,'LineStyle','--')

c=colorbar('southoutside');
% h = colorbar;
% ylabel(h, 'asdasd')
colormap(gca,brewermap([],'*RdBu'))
% colormap('parula');
% colormap default+
% colormasp jet
caxis([0, 0.5]);
c.Ticks = [0:0.1:0.5];
rot = problemparams.rotationmagnitude;
order = problemparams.quadratureorder;
nquad = problemparams.nquadpoints;
type = problemparams.quadraturetype; types = {'tens', 'octa', 'ico'};
conv = problemparams.convolutionflag;
lowrankflag = problemparams.lowrankflag;
whichrank = problemparams.whichrank;

% use  this title in general 
% title(sprintf('r$_{%1.f}$S$_{%i}$, Q$=%s$, $N_q=%i$, Conv $=%i$, LR $=%i$, R $=%i$', ...
%     rot, order, types{type}, nquad, conv, lowrankflag, whichrank), ...
%     'interpreter', 'latex', 'FontSize', FS);
set(gca, 'FontSize', 1.5*FS)

if type ==2
    title(sprintf('r$_{%1.f}$S$_{%i}$, $N_q=%i$', ...
    rot, order,  nquad), ...
    'interpreter', 'latex', 'FontSize', 2.5*FS);
else
    title(sprintf('S$_{%i}$, $N_q=%i$', ...
     order,  nquad), ...
    'interpreter', 'latex', 'FontSize', 2.5*FS);
end

if rot==999
title('Exact solution', 'interpreter', 'latex', 'FontSize', 2.5*FS);
end    


%   a=get(cb); %gets properties of colorbar
%   a.Position %gets the positon and size of the color bar
% set(cb,'Position',[a(1) a(2) 0.70 0.03])

axis equal
axis tight
xlabel('$x$ ','interpreter','latex','FontSize',2*30)
ylabel('$y$ ','interpreter','latex','FontSize',2*30)
% 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1)+ti(1);
% bottom = outerpos(2)+ti(2);
% ax_width = outerpos(3)-ti(1)-ti(3);
% ax_height = outerpos(4)-ti(2)-ti(4);
% ax.Position = [left, bottom, ax_width, ax_height];


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = PP;

name = strcat('/Linesource_',num2str(rot),'_',num2str(nquad));
drawnow
pause(2)
print(strcat(folder, name), '-dpng', '-r0')

%% FIGURE 2: Cross
figure();
plot(linspace(0, 1.5, Ny/2), rhoNumerical(Ny/2, Ny/2+1:end), '-x', 'LineWidth', 2+LW); hold on
d = diag(rhoNumerical); d = d(Ny/2+1:end);
plot(sqrt(2)*linspace(0, 1.5, Ny/2), d, '-x', 'LineWidth', 2+LW);
axis([0, 1.5, 0, 0.7]);
grid on
plot(rhoExact(501:end, 1), rhoExact(501:end, 2), '-', 'LineWidth', 2+LW)
legend('horizontal', 'diagonal', 'exact', 'Location', 'NorthEast','FontSize',2*FS);

% title(sprintf('r$_{%1.f}$S$_{%i}$, Q$=%s$, $N_q=%i$, Conv $=%i$, LR $=%i$, R $=%i$', ...
%     rot, order, types{type}, nquad, conv, lowrankflag, whichrank), ...
%     'interpreter', 'latex', 'FontSize', FS);
% xlabel('$r =\sqrt{x^2+y^2}$','interpreter','latex')
set(gca, 'FontSize', FS)
if type==2
    title(sprintf('r$_{%1.f}$S$_{%i}$, $N_q=%i$', ...
    rot, order,  nquad), ...
    'interpreter', 'latex', 'FontSize', 2.5*FS);
else
    title(sprintf('as-S$_{%i}$, N$_q$=%i', ...
     order,  nquad), ...
    'interpreter', 'latex', 'FontSize', 2.5*FS);
end


% axis equal
% axis tight
% set(gca, 'FontSize',3*FS)
xlabel('$r=\sqrt{x^2+y^2}$','interpreter','latex','FontSize',2*30)
% ylabel('$y$','interpreter','latex','FontSize',30)
ylabel('$\rho$','interpreter','latex','FontSize',2*30)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = PP;
name = strcat('/LinesourceCut_',num2str(rot),'_',num2str(nquad));
drawnow
pause(2)
print(strcat(folder, name), '-dpng', '-r0')

end
