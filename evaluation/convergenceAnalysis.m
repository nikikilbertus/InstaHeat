close all
dim = 2;
dpprefix = 'dp_';
dpvals = 3:12;
rkprefix = 'rk_';
rkvals = {'008', '01', '02', '03', '05', '06'};
path = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';

%load finest dp852 stuff
name = [path dpprefix num2str(max(dpvals)) '.h5'];
dpt = h5read(name, '/time');
dpphi = h5read(name, '/phi');
dppsi = h5read(name, '/psi');
rkphinorms = zeros(1, length(rkvals));
rkpsinorms = zeros(1, length(rkvals));
dtsf = @(i) str2double(['0.' rkvals{i}]);
dts = zeros(1, length(rkvals));
% loading different rk4 dt values
for i = 1:length(rkvals)
    dts(i) = dtsf(i);
    name = [path rkprefix rkvals{i} '.h5'];
    t = h5read(name, '/time');
    if abs(dpt(end) - t(end)) > 1e-15
       error('comparison at different times'); 
    end
    phi = h5read(name, '/phi');
    psi = h5read(name, '/psi');
    rkphinorms(i) = norm(phi(:, end) - dpphi(:, end));
    rkpsinorms(i) = norm(psi(:, end) - dppsi(:, end));
end

slope = logfit(dts,rkphinorms,'loglog');
title(['phi, fix dp853 tol' num2str(10^-max(dpvals)) ', vary dt in rk4, slope=' num2str(slope)]);
xlabel('dt');
ylabel('norm of (rk4 - dp853) at t=250');
shg;
pause;

slope = logfit(dts,rkpsinorms,'loglog');
title(['psi, fix dp853 tol' num2str(10^-max(dpvals)) ', vary dt in rk4, slope=' num2str(slope)]);
xlabel('dt');
ylabel('norm of (rk4 - dp853) at t=250');
shg;
pause;

%load finest rk stuff
name = [path rkprefix rkvals{1} '.h5'];
rkt = h5read(name, '/time');
rkphi = h5read(name, '/phi');
rkpsi = h5read(name, '/psi');

dpphinorms = zeros(1, length(dpvals));
dppsinorms = zeros(1, length(dpvals));
% loading the different dp852 values data
for i = 1:length(dpvals)
    name = [path dpprefix num2str(dpvals(i)) '.h5'];
    t = h5read(name, '/time');
    if abs(rkt(end) - t(end)) > 1e-15
       error('comparison at different times'); 
    end
    phi = h5read(name, '/phi');
    psi = h5read(name, '/psi');
    dpphinorms(i) = norm(phi(:, end) - rkphi(:, end));
    dppsinorms(i) = norm(psi(:, end) - rkpsi(:, end));
end

tols = 10.^-dpvals;

slope = logfit(tols,dpphinorms,'loglog');
title(['phi, fix rk4 dt=' rkvals{1} ', vary tol in dp, slope=' num2str(slope)]);
xlabel('tolerance');
ylabel('norm of (rk4 - dp853) at t=250');
shg;
pause;

slope = logfit(tols,dppsinorms,'loglog');
title(['psi, fix rk4 dt=0.008, vary tol in dp, slope=' num2str(slope)]);
xlabel('tolerance');
ylabel('norm of (rk4 - dp853) at t=250');
shg;
pause;