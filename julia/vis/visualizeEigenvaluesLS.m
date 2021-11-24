function visualizeEigenvalues

% ev = dlmread('/home/thomas/Work/Projects/LowRank/numericalresults/lowrank_vs_fullrank_LS_HighRes/eigenvalues_ls.txt');
% % ev = ev(2:end,:);
% [nt, nev] = size(ev);
% ev = ev./ev(:,1);
% rho = dlmread('/home/thomas/Work/Projects/LowRank/numericalresults/lowrank_vs_fullrank_LS_HighRes/density.txt');
% % rho = rho(2:end,:);
% 
% % errors of the low rank computation, normalized l2
% err = dlmread('/home/thomas/Work/Projects/LowRank/numericalresults/lowrank_vs_fullrank_LS_HighRes/errors.txt');
% 
% save

load matlab

err = err/max(err);


close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'ls_normalized_fast.gif';
spread
for i=1:nt
    subplot(1,2,1)
    %% EV
    hold off
    semilogy(ev(i,:),'-x','LineWidth',3)
    hold on
    xx = [1,10*(1:length(err)-1)];
    semilogy(xx,err,'r-x','LineWidth',3);
    
    
%     med = median(ev(1:i,:));
%     semilogy(med,'b-d');
%     
%     m = mean(ev(1:i,:));
%     semilogy(m,'b-o');
%     
%     gm = geomean(ev(1:i,:));
%     semilogy(gm,'b-v');
    
%     y = myfun(ev(i,:));
%     semilogy(y,'k-x');
%     
%      y = myfun2(ev(i,:));
%     semilogy(y,'k-x');
    
    legend('current (n)EV','(n)final rel. L_2 error of rank  R comp.','Location','SouthWest')
    %,'median (n)EV so far','arith. mean (n)EV so far',...
     %   'geom. mean (n)EV so far','Location','SouthWest');
    
    
    
    title(sprintf('(normalized)EVs and low rank errors for LS at timestep %i.',i))
    xlabel('EV or Rank R ')
%     ylabel('Normalized EV (by largest EV)')
    if i>=100
        axis([1,nev,10^(-5),1.1])
    elseif i>=50
        axis([1,nev,10^(-7),1.1])
    else
        axis([1,nev,10^(-16),1.1])
    end
    grid on
    set(gca,'FontSize',20);
    
    %% Solution
    subplot(1,2,2);
    r = reshape(rho(i,:),204,204);
    imagesc((r));
    colorbar; %caxis([0,0.5]);
    xticks([]);
    yticks([]);
    title(sprintf('Density for LS at timestep %i.',i))
    set(gca,'FontSize',20);
    colormap jet
    drawnow
    
        
    %% gif stuff
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.02,'WriteMode','append');
    end
end





end



function y = myfun(evs)
y = zeros(size(evs));
for i=1:length(evs)-1
   y(i) = sum(evs(i+1:end));
end
y = y/max(y);


end



function y = myfun2(evs)
y = zeros(size(evs));
for i=1:length(evs)-1
   y(i) = mean(evs(i+1:end));
end
y = y/max(y);


end
% 
% 
% function y = myfun3(evs)
% y = zeros(size(evs));
% for i=1:length(evs)-1
%    y(i) = prod(evs(i+1:end));
% end
% y = y/max(y);
% 
% 
% end
% function y = myfun4(evs)
% y = zeros(size(evs));
% for i=1:length(evs)-1
%    y(i) = geommean(evs(i+1:end));
% end
% y = y/max(y);
% 
% 
% end

