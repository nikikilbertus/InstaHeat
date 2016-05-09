%% Laplacian
N = 64; dim = 3;
aa=-pi; bb=pi;
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

%% powerspectrum check bunch davies
bins = 50; L = 10; match = 2;
meff2 = mass^2 - 9 * H(1)^2 / 4;
% n = size(phi,1);
% ps = mkPowerSpectrum(phi,bins,L);
% dps = mkPowerSpectrum(dphi,bins,L);
ps = phips(:,1);
dps = dphips(:,1);

kmax = sqrt(3) * (n/2) * (2*pi/L);
k = linspace(0,kmax,bins);
k2 = k.^2;
dk = 2*pi/L;
kcut = 0.5 * (floor(n/2) + 1) * dk;
kcut2 = kcut^2;
normfac = (sqrt(2 * dk^3 * pi))^(-1);
ps1 = k .* (k2+meff2).^(-.25);
dps1 = k .* (k2+meff2).^(.25);
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
ps2full = ks .* (ks.^2+meff2).^(-.25); % .* exp(-ks.^2 / kcut2);
dps2full = ks .* (ks.^2+meff2).^(.25); % .* exp(-ks.^2 / kcut2);
ps2 = zeros(1, bins); dps2 = zeros(1, bins);
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

loglog(k,ps,k,ps1,k,ps2)
hold on
loglog(k,dps,k,dps1,k,dps2)
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

%% construct karstens IC from fourier
L=1; N=64; mass = 1e-2;

% read fourier karsten
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten/data_64fourier_3.dat';
raw = importdata(name);
k = 2*pi/L;
ks = raw(:,1);
xi = raw(:,2);
yi = raw(:,3);
zi = raw(:,4);
ph = raw(:,5);

% read real space karsten
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
fac = 6.193e-6; % magic number from karste (see mail)
phi = phi * fac + mean(phika(:));
dphi = dphi * fac + mean(dphika(:));
psi = psi * fac + mean(psika(:));
dpsi = dpsi * fac + mean(dpsika(:));

%% Hamiltonian and momentum karsten vs. code
L=1; N=64;
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten/data_64psi_5.dat';
raw = importdata(name);
name = 'comp_karsten_5_simp2';
evaluate3D
phika = reshape(raw(:,4),N,N,N);
dphika = reshape(raw(:,5),N,N,N);
psika = reshape(raw(:,6),N,N,N);
dpsika = reshape(raw(:,7),N,N,N);
rhoka = mkrho(phika,dphika,psika,a(1),mass,L);
N = N(1);
psi = h5read(name,'/psi');
psi = reshape(psi(:,1),N,N,N);
phi = h5read(name,'/phi');
phi = reshape(phi(:,1),N,N,N);
rho = h5read(name,'/rho');
rho = reshape(rho(:,1),N,N,N);
dphi = h5read(name,'/dphi');
dphi = reshape(dphi(:,1),N,N,N);
dpsi = h5read(name,'/dpsi');
dpsi = reshape(dpsi(:,1),N,N,N);
prn = @(f) num2str(max(abs(f(:))));
disp(' '); disp(' ');
disp(['----comparison, mass=' num2str(mass) ', a=' num2str(a(1)) '----'])
[check, t1, t2, t3] = hamiltonianConstraint(psika, dpsika, a(1), rhoka, L);
disp('karstens psi hamiltonian')
disp(['sum: ' prn(check(:))])
disp(['1  : ' prn(t1(:))])
disp(['2  : ' prn(t2(:))])
disp(['3  : ' prn(t3(:))])
[check, t1, t2] = momentumConstraint(psika, dpsika, phika, dphika, rhoka, L);
disp('karstens psi momentum')
disp(['sum: ' num2str(check)])
disp(['1  : ' num2str(t1)])
disp(['2  : ' num2str(t2)])
[check, t1, t2, t3] = hamiltonianConstraint(psi, dpsi, a(1), rho, L);
disp(' ');
disp('code psi hamiltonian')
disp(['sum: ' prn(check(:))])
disp(['1  : ' prn(t1(:))])
disp(['2  : ' prn(t2(:))])
disp(['3  : ' prn(t3(:))])
[check, t1, t2] = momentumConstraint(psi, dpsi, phi, dphi, rho, L);
disp('code psi momentum')
disp(['sum: ' num2str(check)])
disp(['1  : ' num2str(t1)])
disp(['2  : ' num2str(t2)])
disp(' ')
disp(['l2 and l\infty differences in \psi ' num2str(norm(psi(:)-psika(:))) ', ' prn(psi(:)-psika(:))])

%% power spectrum check karsten
L=1; N=64; bins = 50;
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten/data_64psi_5.dat';
raw = importdata(name);
name = 'comp_karsten_5_simp2';
evaluate3D
N = N(1);
phika = reshape(raw(:,4),N,N,N);
dphika = reshape(raw(:,5),N,N,N);
psika = reshape(raw(:,6),N,N,N);
dpsika = reshape(raw(:,7),N,N,N);
rhoka = mkrho(phika, dphika, psika, a(1), mass, L);
rhormska = sqrt(var(rhoka) ./ mean(rhoka).^2);
kmax = sqrt(3) * (N/2) * (2*pi/L);
k = linspace(1,kmax,bins);
phikaps = mkPowerSpectrum(phika,bins,L);
rhormskaps = mkPowerSpectrum(rhormska,bins,L);
loglog(1:bins,phikaps, 1:bins,rhormskaps); shg;

%% Hamiltonian and momentum for bunch davies
nums = [1 2 3 4 5 6 7 8];
L = 10;
mabs = @(f) max(abs(f(:)));
prn = @(f) num2str(max(abs(f(:))));
herrs = zeros(4, length(nums));
mxerrs = zeros(3,length(nums));
myerrs = zeros(3,length(nums));
mzerrs = zeros(3,length(nums));
for i = 1:length(nums)
    name = ['bunch' num2str(nums(i))];
    evaluate3D
    N = N(1);
    psi = h5read(name,'/psi');
    psi = reshape(psi(:,1),N,N,N);
    phi = h5read(name,'/phi');
    phi = reshape(phi(:,1),N,N,N);
    rho = h5read(name,'/rho');
    rho = reshape(rho(:,1),N,N,N);
    dphi = h5read(name,'/dphi');
    dphi = reshape(dphi(:,1),N,N,N);
    dpsi = h5read(name,'/dpsi');
    dpsi = reshape(dpsi(:,1),N,N,N);
    disp(' '); disp(' ');
    disp(['----comparison, mass=' num2str(mass) ', a=' num2str(a(1)) '----'])
    disp('code psi hamiltonian')
    [check, t1, t2, t3] = hamiltonianConstraint(psi, dpsi, a(1), rho, L);
    disp(['sum: ' prn(check(:))])
    disp(['1  : ' prn(t1(:))])
    disp(['2  : ' prn(t2(:))])
    disp(['3  : ' prn(t3(:))])
    herrs(:,i) = [mabs(check) mabs(t1) mabs(t2) mabs(t3)];
    disp('code psi momentum')
    [check, t1, t2] = momentumConstraint(psi, dpsi, phi, dphi, rho, L);
    disp(['sum: ' num2str(check)])
    disp(['1  : ' num2str(t1)])
    disp(['2  : ' num2str(t2)])
    mxerrs(:,i) = [check(1) t1(1) t2(1)];
    myerrs(:,i) = [check(2) t1(2) t2(2)];
    mzerrs(:,i) = [check(3) t1(3) t2(3)];
end
herr = herrs(1,:) ./ max(herrs(2:4,:));
mxerr = mxerrs(1,:) ./ max(mxerrs(2:3,:));
myerr = myerrs(1,:) ./ max(myerrs(2:3,:));
mzerr = mzerrs(1,:) ./ max(mzerrs(2:3,:));
ms = 5 * 10.^(-nums);
loglog(ms, ms.^2 / ms(end)^2 * herrs(1,end),'--', ms, herrs(1,:),'linewidth',2);
xlabel('planck mass'); ylabel('max. abs error of hamiltonian'); shg
figure
loglog(ms, ms.^2 / ms(end)^2 * mxerrs(1,end),'--', ms, mxerrs(1,:),'linewidth',2);
hold on
loglog(ms, myerrs(1,:),'linewidth',2);
loglog(ms, mzerrs(1,:),'linewidth',2);
hold off
xlabel('planck mass'); ylabel('max. abs error of momenutm');

%% masses study for bunch davies
figure
ml = 3; m = 5;
% masses = 1./(2*10.^(-(ml:m)));
masses = 5*10.^(-(ml:m));
rhormsi = zeros(m-ml+1,4);
maxrhos = zeros(m-ml+1,1);
for i = 1:m-ml+1
    name = ['resolutions9/96_5e-' num2str(i+ml-1) '_2e4'];
    evaluate3D
    maxrhos(i) = max(rhorms);
    loglog(a,rhorms/rhorms(1),'linewidth',2); hold on
    legendinfo{i} = ['mass=' num2str(masses(i))];
    rhormsi(i,1) = rhorms(1);
    [~,idx] = min((a - 1e2).^2);
    rhormsi(i,2) = rhorms(idx);
    [~,idx] = min((a - 1e3).^2);
    rhormsi(i,3) = rhorms(idx);
    rhormsi(i,4) = rhorms(end);
end
hold off; xlabel('a'); ylabel('std(\rho) / |<\rho>|');
legend(legendinfo,'location','southeast');
figure
loglog(masses, rhormsi, masses, masses/masses(1) * rhormsi(1,1) * 0.9,'--','linewidth',2);
xlabel('mass'); ylabel('std(\rho) / |<\rho>|');
legend('a=1', 'a=1e2', 'a=1e3', ['a=' num2str(a(end))],'linear reference','location','northwest');

%% extrapolate to nonlinear regime
tmp = 2;
name = ['64_1e5_' num2str(tmp)];
evaluate3D
I = (a>100) & (a<1000);
figure
slope = logfit(a(I),rhorms(I), 'loglog'); shg; pause;
rhormsnonlin = max(maxrhos);
anonlin = a(end) * (rhormsnonlin/rhorms(end))^(1/slope)
aratio = anonlin / a(end);
tnonlin = aratio^(3/2) * t(end)
atmp = linspace(a(find(I,1)), anonlin, 10);
loglog(a, rhorms, atmp, slope * atmp / (slope*atmp(1)) * rhorms(find(I,1))); shg;
estimatedruntime = h5read(name,'/runtime_total') / 3600 * tnonlin/t(end)

%% playing with quantities in karstens paper
L = 10; N = 32;
mpl = 1;
k = 2*pi/L;
alpha = 1;
ms = mass;
Hs = H;
Hend = Hs(1); aend = a(1); Trh = 1e7;
lc = 1./sqrt(3*Hs*ms);
lcfit = fit(log(a), log(mpl*lc)', 'linear');

kphys = k ./ a;
kmingrid = k;
kmaxgrid = sqrt(3) * N/2 * k;
kmaxfit = fit(log(a), log(a/kmaxgrid), 'linear');

% slope = logfit(a, lc, 'loglog'); shg; pause; slope
% slope = logfit(a, 1./H, 'loglog'); shg; pause; slope
% shg; pause;
% ks = [0.02 * a, 0.16 * a, a, 10 * a, 100 * a ];
% h=loglog(a, ks, '--k','linewidth',0.5); hold on;
% arrayfun(@(x) set(get(get(x,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'),h);

kmax = Hend * 1.37e3 * sqrt(ms / (1.4e-6 * mpl)) * (Trh/1e7)^(-1/3) * (Hend/1e13)^(-1/3);
kmin = Hend * 2.74e-6 * (Trh/1e7)^(2/3) * (Hend/1e13)^(-1/3);

deltak = 2/5 * (k^2 ./ (a'.^2 .* Hs.^2) + 3) .* sqrt(psivar);
deltak = deltak / deltak(end) * rhorms(end);

loglog(a, mpl*lc, a, mpl./Hs, a, mpl./kphys, a, rhorms, a, deltak, a, a/kmaxgrid, a, a/kmingrid); hold on; shg;
legend('mpl/sqrt(3 H m)', 'mpl/H', 'mpl/k_{phys}', '\delta \rho / \rho', '\delta k (norm)', 'kmax', 'kmin');
pause;
aint = exp(fminsearch(@(x) (kmaxfit(x) - lcfit(x)).^2, a(end)))
aext = linspace(log(a(1)),log(aint),100);
loglog(exp(aext), exp(kmaxfit(aext)), exp(aext), exp(lcfit(aext)),'linewidth',0.8); hold off; shg;

%% power spectrum analysis
L=10; nbins = 60;
k = 2*pi/L;
name = 'resolutions9/96_5e-3_2e4';
evaluate3D
N = N(1);
lc = 1./sqrt(3*H*mass);
kmin = k;
kmax = sqrt(3) * N/2 * k;
bins = (1:nbins)/nbins * kmax;
for ii = [1 logspace(1,log10(length(phips(1,:))),300)]
    i = int64(floor(ii));
    subplot(2,2,1)
    loglog(bins, phips(1:nbins,i)); xlabel('k'); ylabel('power'); hold on
    k = lc(i)
    kmax = 1/H(i)
    [min(bins) max(bins)]
    lowb = min(phips(1:nbins,i)) + 1e-13;
    higb = max(phips(1:nbins,i));
    [lowb higb]
    loglog([k k],[lowb higb]);
    loglog([kmax kmax],[lowb higb]);
    hold off;
    title(['\phi, a = ' num2str(a(i)) ' / ' num2str(a(end))]);
    subplot(2,2,2)
    plot(a, rhorms, a(i), rhorms(i),'or'); xlabel('a'); ylabel('std \rho / <\rho>');
    shg;
    subplot(2,2,3)
    loglog(bins, rhops(1:nbins,i)); xlabel('k'); ylabel('power'); title('\rho');
    subplot(2,2,4)
    loglog(bins, psips(1:nbins,i)); xlabel('k'); ylabel('power'); title('\psi');
    if i == 1
        Hk = H(i);
        lowb = min(rhops(1:nbins,i));
        higb = max(rhops(1:nbins,i));
        subplot(2,2,3); hold on; plot([Hk Hk], [lowb higb]); hold off;
        lowb = min(psips(1:nbins,i));
        higb = max(psips(1:nbins,i));
        subplot(2,2,4); hold on; plot([Hk Hk], [lowb higb]); hold off;
        pause;
    else
        pause(0.05);
    end
end

%% plot long time bunch davies
name = 'longruns/32_5e-4_1e6_beta_0';
evaluate3D
scal = ones(size(a))';
% scal = a.^(3/2)';
plot(a,phimax.*scal,a,phimin.*scal,a,phimean.*scal); xlabel('a'); ylabel('<\phi>'); shg; pause;
plot(a,dphimax.*scal,a,dphimin.*scal,a,dphimean.*scal); xlabel('a'); ylabel('<d\phi>'); shg; pause;
plot(a,psimax.*scal,a,psimin.*scal,a,psimean.*scal); xlabel('a'); ylabel('<\psi>'); shg; pause;
plot(a,dpsimax.*scal,a,dpsimin.*scal,a,dpsimean.*scal); xlabel('a'); ylabel('<d\psi>'); shg; pause;
plot(a,sqrt(phivar).*scal); xlabel('a'); ylabel('std \phi'); shg; pause;
plot(a,sqrt(dphivar).*scal); xlabel('a'); ylabel('std d\phi'); shg; pause;
plot(a,sqrt(psivar).*scal); xlabel('a'); ylabel('std \psi'); shg; pause;
plot(a,sqrt(dpsivar).*scal); xlabel('a'); ylabel('std d\psi'); shg; pause;

%% long runs
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/longruns/32_5e-4_1e6_beta_004.h5';
a = h5read(name, '/a');
t = h5read(name, '/time');
rhosmry = h5read(name, '/rho_summary');
rhorms = sqrt(rhosmry(2,:) ./ rhosmry(1,:).^2);
figure(1);
semilogx(a,rhorms); xlabel('a'); ylabel('rhorms'); shg;
rhops = h5read(name, '/rho_power_spectrum');
I = (1:1000:length(a));
figure(2);
surf(a(I),1:50,log10(rhops(:,I))); xlabel('a'); ylabel('bins'); shading interp; view(2); title('ps rho');
shg;
figure(3);
plot(diff(t)); hold on
steps = h5read(name, '/steps_total')
a(end)

%% tolerances analysis
base = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/tolerances2/64_5e-3_1e4_tol_';
rtol = 4:2:10; atol = 6:2:16;
relval = 10.^(-rtol); absval = 10.^(-atol);
time = zeros(length(rtol), length(atol));
steps = zeros(size(time));
errinf = zeros(size(time));
phimeanerrl2 = zeros(size(time));
phistderrl2 = zeros(size(time));
as = zeros(size(time));
cstrl2 = zeros(size(time));
getname = @(x,y) [base num2str(x) '_' num2str(y) '.h5'];
phiref = h5read(getname(max(rtol),max(atol)), '/phi_summary');
aref = h5read(getname(max(rtol),max(atol)), '/a');
cstrl2ref = h5read(getname(max(rtol),max(atol)),'/constraints');
cstrl2ref = cstrl2ref(1,end);
% phiref = phiref(1,:);
arefs = h5read(getname(min(rtol),min(atol)), '/a');
phimeanrefs = spline(aref,phiref(1,:),arefs);
phistdrefs = spline(aref,sqrt(phiref(2,:)),arefs);
for i = 1:length(rtol)
    for j = 1:length(atol)
        name = getname(rtol(i),atol(j));
        tols = h5read(name, '/tolerances');
        if relval(i) ~= tols(1) || absval(j) ~= tols(2)
            error('didnt load the right file')
        end
        time(i,j) = h5read(name,'/runtime_stepper');
        steps(i,j) = h5read(name,'/steps_total');
        phi = h5read(name, '/phi_summary');
        a = h5read(name, '/a');
        cstr = h5read(name,'/constraints');
%         semilogy(a,cstr(1,:)); shg; pause;
        cstrl2(i,j) = -log10(abs((cstr(1,end) - cstrl2ref) / cstrl2ref));
        as(i,j) = -log10(abs((a(end) - aref(end))/aref(end)));
%         I = (a>0.9*aref(end));
%         Iref = (aref>0.9*aref(end));
%         plot(a(I),phi(1,I).*a(I)'.^(3/2),aref(Iref),phiref(Iref).*aref(Iref)'.^(3/2)); shg; pause;
        errinf(i,j) = abs((phiref(1,end) - phi(1,end))/phiref(1,end));
        phimean = spline(a, phi(1,:),arefs);
        phistd= spline(a, sqrt(phi(2,:)),arefs);
        phimeanerrl2(i,j) = norm(phimeanrefs - phimean);
        phistderrl2(i,j) = norm(phistdrefs - phistd);
    end
end
% semilogx(relval, time, 'linewidth',2); xlabel('rel tol'); ylabel('time [s]');
% legend('abs: 1e-6', 'abs: 1e-8', 'abs: 1e-10', 'abs: 1e-12'); shg;
% figure
% semilogx(relval, steps, 'linewidth',2); xlabel('rel tol'); ylabel('#steps');
% legend('abs: 1e-6', 'abs: 1e-8', 'abs: 1e-10', 'abs: 1e-12'); shg;
% figure
% semilogx(absval, time', 'linewidth',2); xlabel('abs tol'); ylabel('time [s]');
% legend('rel: 1e-6', 'rel: 1e-8', 'rel: 1e-10', 'rel: 1e-12'); shg;
% figure
% semilogx(absval, steps', 'linewidth',2); xlabel('abs tol'); ylabel('#steps');
% legend('rel: 1e-6', 'rel: 1e-8', 'rel: 1e-10', 'rel: 1e-12'); shg;

figure
bar3(as); set(gca,'XTickLabel',absval); set(gca,'YTickLabel',relval);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 (a_{f} - a_{f}^{ref}) / a_{f}^{ref}');
figure
bar3(steps); set(gca,'XTickLabel',absval); set(gca,'YTickLabel',relval);
xlabel('atol'); ylabel('rtol'); zlabel('#steps');
figure
bar3(cstrl2); set(gca,'XTickLabel',absval); set(gca,'YTickLabel',relval);
xlabel('atol'); ylabel('rtol'); zlabel('hamiltonian constraint norm');
figure
bar3(-log10(errinf)); set(gca,'XTickLabel',absval); set(gca,'YTickLabel',relval);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 (\phi_{f}^{ref} - \phi_{f}) / \phi_{f}^{ref}');
figure
bar3(-log10(phimeanerrl2)); set(gca,'XTickLabel',absval); set(gca,'YTickLabel',relval);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 <\phi> error l_{2}');
figure
bar3(-log10(phistderrl2)); set(gca,'XTickLabel',absval); set(gca,'YTickLabel',relval);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 std \phi error l_{2}');

%% resolutions study
close all
res = [32 48 64 96];
rhos = zeros(length(res),1);
disp('         grid      steps')
for i = 1:length(res)
    name = ['resolutions9/' num2str(res(i)) '_5e-3_2e4'];
    evaluate3D
    disp([N(1) steps])
    legendinfo{i} = num2str(res(i));
    rhos(i) = rhorms(end);
    figure(1)
    plot(diff(t)); hold on
    figure(2)
    loglog(a, rhorms); hold on
    figure(3)
    subplot(4,1,1)
    loglog(a,max(abs(phimin),abs(phimax)),a,abs(phimean)); xlabel('a'); ylabel('<\phi>'); hold on
    subplot(4,1,2)
    loglog(a,sqrt(phivar)); xlabel('a'); ylabel('std \phi'); hold on
    subplot(4,1,3)
    loglog(a,max(abs(dphimin),abs(dphimax)),a,dphimean); xlabel('a'); ylabel('<d\phi>'); hold on
    subplot(4,1,4)
    loglog(a,sqrt(dphivar)); xlabel('a'); ylabel('std d\phi'); hold on
    figure(4)
    subplot(4,1,1)
    loglog(a,psimin,a,psimax); xlabel('a'); ylabel('min/max \psi'); hold on
    subplot(4,1,2)
    loglog(a,sqrt(psivar)); xlabel('a'); ylabel('std \psi'); hold on
    subplot(4,1,3)
    loglog(a,dpsimin,a,dpsimax); xlabel('a'); ylabel('min/max \psi'); hold on
    subplot(4,1,4)
    loglog(a,sqrt(dpsivar)); xlabel('a'); ylabel('std d\psi'); hold on
    figure(5)
    loglog(a, hamcstrl2/res(i)^3); xlabel('a'); ylabel('ham cstr l2'); hold on
    figure(6)
    loglog(a, hamcstrinf); xlabel('a'); ylabel('ham cstr \infty'); hold on
%     figure(7)
%     surf(a,1:50,log10(rhops)); xlabel('a'); ylabel('bins'); shading interp; view(2); title(num2str(res(i)));
%     pause;
end
hold off
figure(1)
xlabel('steps'); ylabel('dt'); legend(legendinfo); shg;
figure(2)
xlabel('a'); ylabel('std \rho / <|\rho|>'); legend(legendinfo); shg;
inflmass
figure
plot(res.^3, rhos, '-o'); xlabel('N^3'); ylabel('final rhorms'); shg;

%% gw power spectrum
close all
num = 50; aup = 500; bins = 60;
name = 'gw1/64_5e-3_2e4';
evaluate3D
gwps = h5read(name, '/gravitational_wave_spectrum');
gwps(gwps < 0) = 0;
imax = find(a>aup,1);
skip = int64(floor(imax/num));
c = 1;
J = jet;
subplot(1,2,2);
loglog(a,rhorms); hold on;
for i = 1:skip:imax
    subplot(1,2,1)
    semilogy(1:bins,gwps(:,i),'color',J(c,:)); hold on
    linfo{c} = ['i = ' num2str(c)]; c=c+1;
    subplot(1,2,2)
    loglog(a(i),rhorms(i),'x','color',J(c,:));
end
hold off; shg;

%% time scaling

tfactor = @(Nold, Nnew) Nnew^3*log(Nnew^3) / (Nold^3*log(Nold^3));