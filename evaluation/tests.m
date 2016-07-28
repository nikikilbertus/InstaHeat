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
L = 10; N = 96;
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
deltak = deltak / deltak(1) * rhorms(1);

loglog(a, mpl*lc, a, mpl./Hs, a, mpl./kphys, a, rhorms, a, deltak, a, a/kmaxgrid, a, a/kmingrid); hold on; shg;
legend('mpl/sqrt(3 H m)', 'mpl/H', 'mpl/k_{phys}', '\delta \rho / \rho', '\delta k (norm)', 'kmax', 'kmin');
% pause;
% aint = exp(fminsearch(@(x) (kmaxfit(x) - lcfit(x)).^2, a(end)))
% aext = linspace(log(a(1)),log(aint),100);
% loglog(exp(aext), exp(kmaxfit(aext)), exp(aext), exp(lcfit(aext)),'linewidth',0.8); hold off; shg;