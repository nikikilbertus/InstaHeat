function [check, t1, t2] = momentumConstraint(psi, dpsi, phi, dphi, H)
N = 64;
dim = 3;
aa=-pi;
bb=pi;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/(bb-aa);
[kx,~,~] = meshgrid(k,k,k);

t1 = ifftn(1i * kx .* fftn(dpsi + H * psi));
phix = ifftn(1i * kx .* fftn(phi));
t2 = - 0.5 .* dphi .* phix;
% psilap = del2(psi, dx) * 2 * dim;
check = t1 + t2;
end