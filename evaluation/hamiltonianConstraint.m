function [check, t1, t2, t3] = hamiltonianConstraint(psi, dpsi, a, rho, L)
N = 64;
dim = 3;
aa=0;
bb=L;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
k = [0:N/2 -N/2+1:-1] * 2*pi/(bb-aa);
[kx,ky,kz] = ndgrid(k,k,k);
kk = kx.^2 + ky.^2 + kz.^2;

H = sqrt(mean(rho(:)) / 3);

psilap = ifftn(-kk.*fftn(psi));
% psilap = del2(psi, dx) * 2 * dim;
t1 = psilap / a^2;
t2 = - 3 * H * (H * psi + dpsi);
t3 = - 0.5 * (rho - mean(rho(:)));
check = t1 + t2 + t3;
end