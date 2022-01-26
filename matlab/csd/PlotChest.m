function PlotChest(x,y,Dose,Rho)
% 2D plotting routine: CT and contour lines.

% Plot density.
h = imagesc(x,y,100*Rho');
set(h,'alphadata',~isnan(Rho'))
xlabel('x [cm]'); ylabel('y [cm]');
axis xy equal tight
colormap(gray)
hold on

% Choose values of contour lines.
Contours = [0.05 0.1 0.25 0.5 0.7 0.8 .9 .95 .98];

% Choose colormap and draw legend
CM = hsv(length(Contours));
for k = 1:length(Contours)
    plot(0,20,'color',CM(k,:))
end
legend(num2str(Contours'*100,'%g%%'));

% Plot isolines.
Dose = Dose/max(max(Dose)); % Normalize dose.
D = Dose';
DoseMin = 0;%min(min(D));
DoseMax = 1;%max(max(D));

for i=1:length(Contours)
    contour(x,y,D,[DoseMin+Contours(i)*(DoseMax-DoseMin) DoseMin+Contours(i)*(DoseMax-DoseMin)],'Color',CM(i,:),'LineWidth',2);
%     h=contour(x,y,100*((z))'+2000,'LineWidth',2);
end

%contour(x,y,D,[-0.03 -0.03],'Color','k','LineWidth',2);

axis xy equal tight
xlabel('x [cm]'), ylabel('y [cm]')