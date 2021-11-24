function visualizeOverview(prefix,plottingparams)



%% plot overview
SigmaS = dlmread(strcat(prefix,'/data/sigmaS.txt'));
SigmaT = dlmread(strcat(prefix,'/data/sigmaT.txt'));
SigmaA = SigmaT-SigmaS;

Source = dlmread(strcat(prefix,'/data/source.txt'));
rho0 = dlmread(strcat(prefix,'/data/rhofinal.txt'));


figure()
subplot(3,3,1); imagesc(SigmaS); colorbar; title('SigmaS');
subplot(3,3,2); imagesc(SigmaA); colorbar; title('SigmaA');
subplot(3,3,3); imagesc(Source); colorbar; title('Source');
subplot(3,3,4); imagesc(rho0); colorbar; title('rho0');


R = dlmread(strcat(prefix,'/data/rhofinal.txt'));
subplot(3,3,5); imagesc(R); colorbar; colormap jet; title('final density');
subplot(3,3,6); imagesc(log10(abs(R)));colorbar; caxis([-7,0]); colormap jet;  title('final density, log10');


p = dlmread(strcat(prefix,'/data/quadpoints.txt'));
w = dlmread(strcat(prefix,'/data/quadweights.txt'));

subplot(3,3,7); plot3(p(:,1),p(:,2),p(:,3),'x'); title('quadpoints, xyz');
subplot(3,3,8); plot(p(:,1),p(:,2),'x'); title('quadpoints, xy');
subplot(3,3,9); histogram(w,10); title('weight distribution');



fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = plottingparams.PP;
print(strcat(prefix,'/Overview'), '-dpng', '-r0')


end
