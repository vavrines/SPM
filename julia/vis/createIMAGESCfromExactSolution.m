function createIMAGESCfromExactSolution


data = dlmread('exactLineSource.txt');
r = data(:,1);
rho = data(:,2);


x = linspace(-1.5,1.5,204); y = x;

[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);

Rlin = R(:);
Philin = interp1(r,rho,Rlin);
Phi = reshape(Philin,size(R));
Phi(isnan(Phi)) = 0.0;

close all
figure()
imagesc(Phi)
dlmwrite('exactRho.txt',Phi);










end