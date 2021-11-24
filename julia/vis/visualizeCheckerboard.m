function rhoNumerical = visualizeCheckerboard(folder, plottingparams, problemparams)
close all
FS = plottingparams.FS; % font size
LW = plottingparams.LW; % line width
PP = plottingparams.PP; % paper position

rhoNumerical = flipud(dlmread(strcat(folder, '/data/rhofinal.txt')));

rhoNumerical = rhoNumerical(3:end-2, 3:end-2); % remove ghost cells

[Ny, Nx] = size(rhoNumerical);
X = linspace(0, 7, Nx); Y = linspace(0, 7, Ny);

%% Figure 1: IMAGESC
figure()
rhoNumerical(rhoNumerical < 0) = 1e-10;
h = pcolor(X, Y, log10(abs(flipud(rhoNumerical))));
set(h, 'EdgeColor', 'none');
hold on
line([1, 1], [0, 7], 'Color', 'red', 'LineWidth', 2*LW); %,'LineStyle','--')
line([0, 7], [6, 6], 'Color', 'blue', 'LineWidth', 2*LW); %,'LineStyle','--')
xlim([0, 7]);
ylim([0, 7]);
xticks([]);
yticks([]);
c=colorbar('southoutside')
colormap('jet');
caxis([-7, 0]);
c.Ticks = [-7:1:0];
rot = problemparams.rotationmagnitude;
order = problemparams.quadratureorder;
nquad = problemparams.nquadpoints;
type = problemparams.quadraturetype; types = {'tens', 'octa', 'ico'};
conv = problemparams.convolutionflag;
lowrankflag = problemparams.lowrankflag;
whichrank = problemparams.whichrank;
% title(sprintf('r$_{%1.f}$S$_{%i}$, Q$=%s$, $N_q=%i$, Conv $=%i$, LR $=%i$, R $=%i$', ...
%     rot, order, types{type}, nquad, conv, lowrankflag, whichrank), ...
%     'interpreter', 'latex', 'FontSize', FS);
set(gca, 'FontSize', 1.5*30)

if type==2
    title(sprintf('r$_{%1.f}$S$_{%i}$, $N_q=%i$', ...
    rot, order,  nquad), ...
    'interpreter', 'latex', 'FontSize', 2.5*FS);
else
    title(sprintf('S$_{%i}$, $N_q=%i$', ...
     order,  nquad), ...
    'interpreter', 'latex', 'FontSize', 2.5*FS);
end

axis equal
axis tight
hold on
xlabel('$x$','interpreter','latex','FontSize',2*30)
ylabel('$y$','interpreter','latex','FontSize',2*30)

% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1)+ti(1);
% bottom = outerpos(2)+ti(2);
% ax_width = outerpos(3)-ti(1)-ti(3);
% ax_height = outerpos(4)-ti(2)-ti(4);
% ax.Position = [left, bottom, ax_width, ax_height];
% drawnow
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = PP;
name = strcat('/Checkerboard',num2str(rot),'_',num2str(nquad));
drawnow
pause(2)
print(strcat(folder, name), '-dpng', '-r0')

% %% FIGURE 2: Cross
% figure();
% 
% F = griddedInterpolant({X, Y}, rhoNumerical, 'cubic');
% yq = Y; xq = 1;
% zhori = F({xq, yq});
% plot(yq, log10(zhori), '-x', 'LineWidth', LW); hold on
% yq = 1; xq = X;
% zverti = F({xq, yq});
% grid on
% plot(xq, log10(zverti), '-x', 'LineWidth', LW); hold on
% legend('horizontal', 'vertical', 'Location', 'SouthWest','FontSize',2*FS);
% % title(sprintf('r$_{%1.f}$S$_{%i}$, Q$=%s$, $N_q=%i$, Conv $=%i$, LR $=%i$, R $=%i$', ...
% %     rot, order, types{type}, nquad, conv, lowrankflag, whichrank), ...
% %     'interpreter', 'latex', 'FontSize', FS);
% set(gca, 'FontSize', FS)
% if type==2
%     title(sprintf('r$_{%1.f}$S$_{%i}$, $N_q=%i$', ...
%     rot, order,  nquad), ...
%     'interpreter', 'latex', 'FontSize', 2.5*FS);
% else
%     title(sprintf('S$_{%i}$, $N_q=%i$', ...
%      order,  nquad), ...
%     'interpreter', 'latex', 'FontSize', 2.5*FS);
% end
% xlabel('$x$ (vertical) or $y$ (horizontal)','interpreter','latex','FontSize',2*30)
% % ylabel('$y$','interpreter','latex','FontSize',30)
% ylabel('log${}_{10}(\rho)$','interpreter','latex','FontSize',2*30)
% 
% axis([0, 7, -7, 0])
% % set(gca, 'FontSize', FS)
% 
% 
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = PP;
% name = strcat('/CheckerboardCut',num2str(rot),'_',num2str(nquad));
% drawnow
% pause(2)
% print(strcat(folder, name), '-dpng', '-r0')

end
