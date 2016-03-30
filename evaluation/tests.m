% Laplacian

N = 64;
dim = 3;
aa=-pi;
bb=pi;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
[xx,yy,zz] = meshgrid(x,x,x);

% test = exp(-(xx.^2 + yy.^2 + zz.^2));
% testl = (4*(xx.^2+yy.^2+zz.^2)-6).*test;
test = sin(xx) .* sin(yy) .* sin(zz);
testl = -3 * sin(xx) .* sin(yy) .* sin(zz);
testg = cos(xx).^2 .* sin(yy).^2 .* sin(zz).^2 + ... 
        sin(xx).^2 .* cos(yy).^2 .* sin(zz).^2 + ...
        sin(xx).^2 .* sin(yy).^2 .* cos(zz).^2;

k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/(bb-aa);
k2 = [0:N/2 -N/2+1:-1] * 2*pi/(bb-aa);
[X,Y,Z] = meshgrid(k2,k2,k2);
kk = X.^2 + Y.^2 + Z.^2;
[kx,ky,kz] = meshgrid(k,k,k);
phix = ifftn(1i * kx .* fftn(test));
phiy = ifftn(1i * ky .* fftn(test));
phiz = ifftn(1i * kz .* fftn(test));

testg1 = phix.^2 + phiy.^2 + phiz.^2;
testl1 = ifftn(-kk.*fftn(test));
testl2 = del2(test, dx) * 2 * dim;