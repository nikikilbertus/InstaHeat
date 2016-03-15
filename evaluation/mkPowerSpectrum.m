function [ps] = mkPowerSpectrum(f, bins)

N = numel(f);
fk = fftn(f);
fk2 = abs(fk).^2;

if ndims(f) == 1
    ks = sqrt([0:N/2 -N/2+1:-1].^2);
elseif ndims(f) == 2
    [Nx, Ny] = size(f);
    kx = [0:Nx/2 -Nx/2+1:-1];
    ky = [0:Ny/2 -Ny/2+1:-1];
    [X,Y] = meshgrid(kx,ky);
    ks = sqrt(X.^2 + Y.^2);
elseif ndims(f) == 3
    [Nx, Ny, Nz] = size(f);
    kx = [0:Nx/2 -Nx/2+1:-1];
    ky = [0:Ny/2 -Ny/2+1:-1];
    kz = [0:Nz/2 -Nz/2+1:-1];
    [X,Y,Z] = meshgrid(kx,ky,kz);
    ks = sqrt(X.^2 + Y.^2 + Z.^2);
end

kmax = max(ks(:));
ps = zeros(1, bins);
for i = 1:N
    idx = int64(fix(bins * (ks(i) / kmax) - 1e-10) + 1);
    if(idx > bins)
        error('wrong index');
    end
    ps(idx) = ps(idx) + fk2(i)/N;
end
% ps0 = fk2(1)/N^2;

end

%%%%%%%%%%%%%%%%%%%%
% [~, I] = sort(kks(:));
% fmom = fmom(I);
% perbin = ceil(Nall / bins)
% newperbin = perbin;
% ps = zeros(1, bins);
% inds = zeros(1,Nall);
% for i = 1:bins
%     tmp = 0;
%     if i == bins && bins*perbin > Nall
%         newperbin = bins * perbin - Nall;
%         disp('perbin set to')
%         disp(newperbin)
%     end
%     for j = 1:min(perbin, newperbin)
%         tmp = tmp + fmom((i-1)*perbin + j);
%         inds((i-1)*perbin + j) = (i-1)*perbin + j;
%     end
%     ps(i) = tmp;
% end