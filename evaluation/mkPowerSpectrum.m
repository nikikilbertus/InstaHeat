function [ps, ps0] = mkPowerSpectrum(f, bins, L)
    N = numel(f);
    fk = fftn(f);
    fk2 = abs(fk).^2;

    if ndims(f) == 1
        ks = sqrt([0:N/2 -N/2+1:-1].^2) * 2*pi/L;
    elseif ndims(f) == 2
        [Nx, Ny] = size(f);
        kx = [0:Nx/2 -Nx/2+1:-1];
        ky = [0:Ny/2 -Ny/2+1:-1];
        [X,Y] = meshgrid(kx,ky);
        ks = sqrt(X.^2 + Y.^2) * 2*pi/L;
    elseif ndims(f) == 3
        [Nx, Ny, Nz] = size(f);
        kx = [0:Nx/2 -Nx/2+1:-1];
        ky = [0:Ny/2 -Ny/2+1:-1];
        kz = [0:Nz/2 -Nz/2+1:-1];
        [X,Y,Z] = meshgrid(kx,ky,kz);
        ks = sqrt(X.^2 + Y.^2 + Z.^2) * 2*pi/L;
    end

    kmax = max(ks(:));
    ps = zeros(1, bins);
    % treat 0 mode (i.e. the constant offset) separately
    ps0 = fk2(1)/N;
    for i = 2:N
        idx = int64(fix(bins * (ks(i) / kmax) - 1e-14) + 1);
        if(idx > bins)
            error('the index went out of range, this should not happen!');
        end
        ps(idx) = ps(idx) + fk2(i)/N;
    end
end