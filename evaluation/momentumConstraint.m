function [c, tmp1, tmp2] = momentumConstraint(psi, dpsi, phi, dphi, rho, L)
N = size(psi,1);
dim = ndims(psi);
aa=0;
bb=L;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/(bb-aa);
[kx,ky,kz] = meshgrid(k,k,k);

H = sqrt(mean(rho(:)) / 3);

dir = kx;
t1 = ifftn(1i * dir .* fftn(dpsi + H * psi));
phii = ifftn(1i * dir .* fftn(phi));
t2 = - 0.5 .* dphi .* phii;
% psilap = del2(psi, dx) * 2 * dim;
check = t1 + t2;
c(1) = max(abs(check(:)));
tmp1(1) = max(abs(t1(:)));
tmp2(1) = max(abs(t2(:)));

dir = ky;
t1 = ifftn(1i * dir .* fftn(dpsi + H * psi));
phii = ifftn(1i * dir .* fftn(phi));
t2 = - 0.5 .* dphi .* phii;
% psilap = del2(psi, dx) * 2 * dim;
check = t1 + t2;
c(2) = max(abs(check(:)));
tmp1(2) = max(abs(t1(:)));
tmp2(2) = max(abs(t2(:)));

dir = kz;
t1 = ifftn(1i * dir .* fftn(dpsi + H * psi));
phii = ifftn(1i * dir .* fftn(phi));
t2 = - 0.5 .* dphi .* phii;
% psilap = del2(psi, dx) * 2 * dim;
check = t1 + t2;
c(3) = max(abs(check(:)));
tmp1(3) = max(abs(t1(:)));
tmp2(3) = max(abs(t2(:)));

% modified
% t1 = dpsi + H * psi;
% t2 = - 0.5 * mean(dphi(:)) .* (phi - mean(phi(:)));
% check = t1 + t2;
end