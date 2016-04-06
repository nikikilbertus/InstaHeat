%% Laplacian

N = 64;
dim = 3;
aa=-pi;
bb=pi;
dx = (bb-aa)/N;
x = aa + (0:N-1)*dx;
[xx,yy,zz] = meshgrid(x,x,x);

% test = exp(-(xx.^2 + yy.^2 + zz.^2));
% testl = (4*(xx.^2+yy.^2+zz.^2)-6).*test;
test = sin(xx) .* sin(yy) .* sin(zz);
testl = -3 * sin(xx) .* sin(yy) .* sin(zz);
testg = cos(xx).^2 .* sin(yy).^2 .* sin(zz).^2 + ... 
        sin(xx).^2 .* cos(yy).^2 .* sin(zz).^2 + ...
        sin(xx).^2 .* sin(yy).^2 .* cos(zz).^2;

k = [0:N/2-1 0 -N/2+1:-1] * 2*pi/(bb-aa);
k2 = [0:N/2 -N/2+1:-1] * 2*pi/(bb-aa);
[X,Y,Z] = meshgrid(k2,k2,k2);
kk = X.^2 + Y.^2 + Z.^2;
[kx,ky,kz] = meshgrid(k,k,k);
phix = ifftn(1i * kx .* fftn(test));
phiy = ifftn(1i * ky .* fftn(test));
phiz = ifftn(1i * kz .* fftn(test));

testg1 = phix.^2 + phiy.^2 + phiz.^2;
testl1 = ifftn(-kk.*fftn(test));
testl2 = del2(test, dx) * 2 * dim;

%% powerspectrum check
bins = 80;
L=10;
n = size(phi,1);
meff2 = mass^2 - 9 * H(1)^2 / 4;
match = 2;

ps = mkPowerSpectrum(phi,bins,L);
dps = mkPowerSpectrum(dphi,bins,L);

kmax = sqrt(3) * (n/2) * (2*pi/L);
k = linspace(0,kmax,bins);
k2 = k.^2;
dk = 2*pi/L;
kcut = 0.5 * (floor(n/2) + 1) * dk;
kcut2 = kcut^2;
normfac = (sqrt(2 * dk^3 * pi))^(-1);
ps1 = k * dk .* (k2+meff2).^(-.25);
dps1 = k * dk .* (k2+meff2).^(.25);
% ps1 = ps1 * normfac .* exp(-k2 / kcut2);
% dps1 = dps1 * normfac .* exp(-k2 / kcut2);
ps1 = ps1 / ps1(match) * ps(match);
dps1 = dps1 / dps1(match) * dps(match);

% analytic spectra with binning
[Nx, Ny, Nz] = size(phi);
kx = [0:Nx/2 -Nx/2+1:-1];
ky = [0:Ny/2 -Ny/2+1:-1];
kz = [0:Nz/2 -Nz/2+1:-1];
[X,Y,Z] = meshgrid(kx,ky,kz);
ks = sqrt(X.^2 + Y.^2 + Z.^2) * 2*pi/L;
ps2full = ks.^(-1) .* (ks.^2+meff2).^(-.25) * normfac .* exp(-ks / kcut2);
dps2full = ks.^(-1) .* (ks.^2+meff2).^(.25) * normfac .* exp(-ks / kcut2);
ps2 = zeros(1, bins);
dps2 = zeros(1, bins);
for i = 2:N
    idx = int64(fix(bins * (ks(i) / kmax) - 1e-10) + 1);
    if(idx > bins)
        error('wrong index');
    end
    ps2(idx) = ps2(idx) + ps2full(i);
    dps2(idx) = dps2(idx) + dps2full(i);
end
ps2 = ps2 / ps2(match) * ps(match);
dps2 = dps2 / dps2(match) * dps(match);

loglog(k,ps,k,ps1)%,k,ps2)
hold on
loglog(k,dps,k,dps1)%,k,dps2)
hold off
shg

%% power spectrum check: total power
bins = 40;
L=10;
n = size(phi,1);
meff2 = mass^2 - 9 * H(1)^2 / 4;
kmax = sqrt(3) * (n/2) * (2*pi/L);
dk = 2*pi/L;
kcut = 0.5 * (floor(n/2) + 1) * dk;
normfac = (sqrt(2 * dk^3 * pi) * n^3)^(-1);
mrat = 5e-6;

[ps,ps0] = mkPowerSpectrum(phi, bins, L);
n1 = sum(ps) + ps0;
n2 = norm(phi(:))^2;

xiint = @(k, r) k.*(k.^2+meff2).^(-.25)*normfac.*exp(-k.^2/kcut2).*sin(k.*r)*mrat./r;
xi = @(r) integral(@(k) xiint(k,r), 0, kmax);


dx = L/n;
x = (0:n-1)*dx;

[X,Y,Z] = ndgrid(x,x,x);
R = sqrt(X.^2 + Y.^2 + Z.^2);
phigen = zeros(size(phi));
for i = 1:numel(phigen)
    phigen(i) = xi(R(i));
end


%% how far is rho0 from rho
tmp = abs( rho0 - mean(rho0(:)) );
mean(tmp(:))
