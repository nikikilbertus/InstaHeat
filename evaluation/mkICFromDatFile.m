%% construct initial conditions from fourier modes in dat file
% construct real space initial conditions from Karstens initial conditions in
% Fourier space, check whether they match

% specified in the dat files
L = 1; N = 64; mass = 1e-2;

% read the Fourier modes
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten/data_64fourier_3.dat';
raw = importdata(name);
k = 2*pi/L;
ks = raw(:,1);
xi = raw(:,2);
yi = raw(:,3);
zi = raw(:,4);
ph = raw(:,5);

% read corresponding real space data
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten/data_64psi_3.dat';
raw2 = importdata(name);
phika = reshape(raw2(:,4),N,N,N);
dphika = reshape(raw2(:,5),N,N,N);
psika = reshape(raw2(:,6),N,N,N);
dpsika = reshape(raw2(:,7),N,N,N);

nn = 16;
indices = (0:16);
[kx,ky,kz] = meshgrid(indices,indices,indices);
ks1 = sqrt(kx.^2 + ky.^2 + kz.^2) * k;
ks2 = zeros(nn,nn,nn);
for i = 1:length(ks)
    ii = int64(xi(i));
    jj = int64(yi(i));
    kk = int64(zi(i));
    ks2(ii+1,jj+1,kk+1) = ks(i);
end

cphi = raw(:,6);
sphi = raw(:,7);
cdphi = raw(:,8);
sdphi = raw(:,9);
cpsi = raw(:,10);
spsi = raw(:,11);
cdpsi = raw(:,12);
sdpsi = raw(:,13);

x = linspace(0,L-L/N,N); y = x; z = x;
[x,y,z] = ndgrid(x,y,z);

phi = zeros(N,N,N);
dphi = zeros(N,N,N);
psi = zeros(N,N,N);
dpsi = zeros(N,N,N);
for i = 1:length(ks)
    in = k * (xi(i) * x + yi(i) * y + zi(i) * z) + ph(i);
    ct = cos(in);
    st = sin(in);
    phi = phi + cphi(i) *  ct - sphi(i) * st;
    dphi = dphi + cdphi(i) *  ct - sdphi(i) * st;
    psi = psi + cpsi(i) *  ct - spsi(i) * st;
    dpsi = dpsi + cdpsi(i) *  ct - sdpsi(i) * st;
    disp(i);
end
fac = 6.193e-6; % magic number as overall factor (from Karsten, see mail)
phi = phi * fac + mean(phika(:));
dphi = dphi * fac + mean(dphika(:));
psi = psi * fac + mean(psika(:));
dpsi = dpsi * fac + mean(dpsika(:));