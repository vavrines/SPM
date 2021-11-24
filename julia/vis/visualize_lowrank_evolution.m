function visualize_lowrank_evolution


% rho = dlmread('/home/thomas/Work/Projects/LowRank/numericalresults/lowrank_vs_fullrank_CB_HighRes/density.txt');
% save('rho','rho');

% d = dir('../out/20*')
% for i = 1:length(d)
%     folder = d(i).name
%     prefix = strcat('../out/', folder);
%     filename = strcat(prefix,'/rho_timeevolution');
% %     LR{i} 
%     rho = dlmread(filename);
%     i
%     save(strcat('rho_',num2str(i)),'rho');
% 
% end
% 
% % save('LR','LR');
% % save('rho','rho');


keyboard

figure()
spread
[nt,nx2] = size(rho);
rho1 = LR{3};

for i=1:nt
    subplot(1,2,1)
    imagesc(log10(reshape(rho(i,:),[284,284])));
    colorbar;  caxis([-7,0]);
    subplot(1,2,2)
    imagesc(log10(reshape(rho1(i,:),[284,284])));
    colorbar; caxis([-7,0]);
    drawnow
end



end
