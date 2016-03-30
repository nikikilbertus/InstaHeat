function [rho] = mkrho(phi, dphi, psi, a, mass)
N = 64;
dim = 3;
aa=-pi;
bb=pi;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/(bb-aa);
[kx,ky,kz] = meshgrid(k,k,k);

phix = ifftn(1i * kx .* fftn(phi));
phiy = ifftn(1i * ky .* fftn(phi));
phiz = ifftn(1i * kz .* fftn(phi));
grad = phix.^2 + phiy.^2 + phiz.^2;
rho = 0.5 * (1 - 2 * psi) .* dphi.^2 + ...
      0.5 * (1 + 2 * psi) .* grad / a^2 + ...
      0.5 * mass^2 * phi.^2;
end