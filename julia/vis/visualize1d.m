function rho = visualize1d(folder, plottingparams, problemparams)

FS = plottingparams.FS; % font size
LW = plottingparams.LW; % line width
PP = plottingparams.PP; % paper position

rhoExact = dlmread('exactLineSource.txt');
rhoNumerical = dlmread(strcat(folder, '/rhofinal.txt'));

rhoNumerical = rhoNumerical(3:end-2, 3:end-2); % remove ghost cells
rho = rhoNumerical(2, :);

%% FIGURE 1: Cross
figure();
plot(rhoNumerical(2, :), '-x', 'LineWidth', LW); hold on
grid on
rot = problemparams.rotationmagnitude;
order = problemparams.quadratureorder;
nquad = problemparams.nquadpoints;
type = problemparams.quadraturetype; types = {'tens', 'octa', 'ico','dim2'};
conv = problemparams.convolutionflag;
lowrankflag = problemparams.lowrankflag;
whichrank = problemparams.whichrank;
title(sprintf('r$_{%1.f}$S$_{%i}$, Q$=%s$, $N_q=%i$, Conv $=%i$, LR $=%i$, R $=%i$', ...
    rot, order, types{type}, nquad, conv, lowrankflag, whichrank), ...
    'interpreter', 'latex', 'FontSize', FS);
set(gca, 'FontSize', FS)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = PP;
print(strcat(folder, '/1d'), '-dpng', '-r0')

end
