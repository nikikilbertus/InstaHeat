% Laplacian
N = 64;
aa=0;
bb=10;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
[xx,yy,zz] = meshgrid(x,x,x);
k = [0:N/2 -N/2+1:-1] * 2*pi/(bb-aa);
[X,Y,Z] = meshgrid(k,k,k);
kk = X.^2 + Y.^2 + Z.^2;
test = exp(-(xx.^2 + yy.^2 + zz.^2));
testl1 = (4*(xx.^2+yy.^2+zz.^2)-6).*test;
testl2 = ifftn(-kk.*fftn(test));
check = testl1 - testl2;