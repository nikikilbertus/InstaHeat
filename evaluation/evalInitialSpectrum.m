%% evaluate power spectrum on initial time slice

%% setup (user input)
name = 'test40';

%% run
require('phi','dphi','phips','N','spatial_bounds_x','rhoS','inflaton_mass','mass');
readDsets;
N = N(1);
nbins = size(phips,2);
L = spatial_bounds_x(2)-spatial_bounds_x(1);
kmin = 2*pi/L; kmax = sqrt(3)*kmin*N/2;
k = (1:nbins)*kmax/nbins;
k2 = k.^2;
meff2 = mass^2-9*H(1)^2/4;
phi = squeeze(phi(1,:,:,:));
dphi = squeeze(dphi(1,:,:,:));
% ps = mkPowerSpectrum(phi,nbins,L);
% dps = mkPowerSpectrum(dphi,nbins,L);
ps = phips(1,:);
% dps = dphips(1,:);

% kcut = 0.5 * (floor(n/2) + 1) * dk;
% kcut2 = kcut^2;
normfac = inflaton_mass / (N^3 * sqrt(2 * pi * kmin^3));
ps1 = normfac * k .* (k2+meff2).^(-.25);
dps1 = normfac * k .* (k2+meff2).^(.25);
% ps1 = ps1 .* exp(-k2 / kcut2);
% dps1 = dps1 .* exp(-k2 / kcut2);
% ps1 = ps1 / ps1(match) * ps(match);
% dps1 = dps1 / dps1(match) * dps(match);

loglog(k,ps,k,ps1); hold on; shg;
% figure
% loglog(k,dps,k,dps1); shg;
return

% analytic spectra with binning
[Nx, Ny, Nz] = size(phi);
kx = [0:Nx/2 -Nx/2+1:-1];
ky = [0:Ny/2 -Ny/2+1:-1];
kz = [0:Nz/2 -Nz/2+1:-1];
[X,Y,Z] = meshgrid(kx,ky,kz);
ks = sqrt(X.^2+Y.^2+Z.^2)*kmin;
ps2full = ks.*(ks.^2+meff2).^(-.25); % .* exp(-ks.^2 / kcut2);
dps2full = ks.*(ks.^2+meff2).^(.25); % .* exp(-ks.^2 / kcut2);
ps2 = zeros(1,bins); dps2 = zeros(1,bins);
for i = 2:N
    idx = int64(fix(nbins*(ks(i)/kmax)-1e-14)+1);
    if(idx > bins)
        error('wrong index');
    end
    ps2(idx) = ps2(idx)+ps2full(i);
    dps2(idx) = dps2(idx)+dps2full(i);
end
ps2 = ps2/ps2(match)*ps(match);
dps2 = dps2/dps2(match)*dps(match);

loglog(k,ps,k,ps1,k,ps2)
hold on
loglog(k,dps,k,dps1,k,dps2)
hold off
shg