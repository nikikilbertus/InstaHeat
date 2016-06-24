function [c, tmp1, tmp2] = mkMomentumConstraints(psi, dpsi, phi, dphi, rho, L)
    if ndims(phi) ~= 3
       error('''mkMomentumConstraints'' only works in three dimensions'); 
    end
    N = size(psi,1);
    k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/L;
    [kx,ky,kz] = meshgrid(k,k,k);
    H = sqrt(mean(rho(:)) / 3);
    dir = kx;
    t1 = ifftn(1i * dir .* fftn(dpsi + H * psi));
    phii = ifftn(1i * dir .* fftn(phi));
    t2 = - 0.5 .* dphi .* phii;
    % t2 = - 0.5 .* (1 - 2 * psi) .* dphi .* phii;
    % psilap = del2(psi, dx) * 2 * dim;
    check = t1 + t2;
    c(1) = max(abs(check(:)));
    tmp1(1) = max(abs(t1(:)));
    tmp2(1) = max(abs(t2(:)));

    dir = ky;
    t1 = ifftn(1i * dir .* fftn(dpsi + H * psi));
    phii = ifftn(1i * dir .* fftn(phi));
    t2 = - 0.5 .* dphi .* phii;
    % t2 = - 0.5 .* (1 - 2 * psi) .* dphi .* phii;
    % psilap = del2(psi, dx) * 2 * dim;
    check = t1 + t2;
    c(2) = max(abs(check(:)));
    tmp1(2) = max(abs(t1(:)));
    tmp2(2) = max(abs(t2(:)));

    dir = kz;
    t1 = ifftn(1i * dir .* fftn(dpsi + H * psi));
    phii = ifftn(1i * dir .* fftn(phi));
    t2 = - 0.5 .* dphi .* phii;
    % t2 = - 0.5 .* (1 - 2 * psi) .* dphi .* phii;
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