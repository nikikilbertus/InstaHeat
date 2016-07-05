function [check, t1, t2, t3] = mkHamiltonianConstraint(psi, dpsi, a, rho, L)
    psi = squeeze(psi);
    dpsi = squeeze(dpsi);
    rho = squeeze(rho);
    if ndims(psi) ~= 3
       error('''mkHamiltonianConstraint'' only works in three dimensions'); 
    end
    N = size(psi,1);
    k = [0:N/2 -N/2+1:-1] * 2*pi/L;
    [kx,ky,kz] = ndgrid(k,k,k);
    kk = kx.^2 + ky.^2 + kz.^2;
    H = sqrt(mean(rho(:)) / 3);
    psilap = ifftn(-kk.*fftn(psi));
    % psilap = del2(psi, dx) * 2 * dim;
    t1 = psilap;
    t2 = - 3 * H * (H * psi + dpsi) .* a^2;
    t3 = - 0.5 * (rho - mean(rho(:))) .* a^2;
    check = t1 + t2 + t3;
end