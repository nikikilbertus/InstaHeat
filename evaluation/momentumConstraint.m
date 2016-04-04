function [check, t1, t2] = momentumConstraint(psi, dpsi, phi, dphi, rho)
N = 64;
dim = 3;
aa=0;
bb=1;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/(bb-aa);
[kx,ky,kz] = meshgrid(k,k,k);

H = sqrt(mean(rho(:)) / 3 );

t1 = ifftn(1i * kx .* fftn(dpsi + H * psi));
phii = ifftn(1i * kx .* fftn(phi));
t2 = - 0.5 .* dphi .* phii;
% psilap = del2(psi, dx) * 2 * dim;
check = t1 + t2;

% modified
% t1 = dpsi + H * psi;
% t2 = - 0.5 * mean(dphi(:)) .* (phi - mean(phi(:)));
% check = t1 + t2;
end