function [check, t1, t2, t3] = hamiltonianConstraint(psi, dpsi, a, H, rho)
N = 64;
dim = 3;
aa=-pi;
bb=pi;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
k = [0:N/2 -N/2+1:-1] * 2*pi/(bb-aa);
[X,Y,Z] = meshgrid(k,k,k);
kk = X.^2 + Y.^2 + Z.^2;

psilap = ifftn(-kk.*fftn(psi));
% psilap = del2(psi, dx) * 2 * dim;
t1 = psilap / a^2;
t2 = - 3 * H * (H * psi + dpsi);
t3 = - 0.5 * (rho - mean(rho(:)));
check = t1 + t2 + t3;
end