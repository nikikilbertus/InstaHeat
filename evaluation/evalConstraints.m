%% DEPRECATED! (constraints are now computed during the simulation)
%evaluate the Hamiltonian and Momentum constraints

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