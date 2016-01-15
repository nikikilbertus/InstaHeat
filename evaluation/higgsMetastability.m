close all
dim = 2;
name = 'dp_12';

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
t = h5read(name, '/time');
a = h5read(name, '/a');
rho = h5read(name, '/rho');
phi = h5read(name, '/phi');
psi = h5read(name, '/psi');
powspec = h5read(name, '/power_spectrum');
% compute some further properties
H = sqrt(rho / 3);
phiAvg = mean(phi);
psiAvg = mean(psi);
Nt = length(t);

% start with some simple plots (rho over a with a^-4 for reference)
% loglog(a, rho);
loglog(a, rho, a, a.^(-4) * rho(1));
% 'logfit' works only with logfit package at
% http://www.mathworks.com/matlabcentral/fileexchange/29545
% (just put the logfit.m file in the same folder as this one)
% logfit(a, rho, 'loglog')
xlabel('a')
ylabel('rho');
legend('data', 'reference: a^{-4}');
shg;
pause;

% H over a with a^-2 for reference
% loglog(a, H);
loglog(a, H, a, a.^(-2) * H(1));
% logfit(a, H, 'loglog')
xlabel('a')
ylabel('H');
legend('data', 'reference: a^{-2}');
shg;
pause;

% average field value and perturbation (with max and min) over time
subplot(1,2,1)
plot(t, phiAvg, t, max(phi), t, min(phi));
title('phi')
xlabel('t')
legend('<phi>', 'max', 'min')
subplot(1,2,2)
plot(t, psiAvg, t, max(psi), t, min(psi));
title('psi')
xlabel('t')
legend('<psi>', 'max', 'min')
shg;
pause;
close all;

% the power spectrum with a reference plane at 10^-10 (everything below
% that might as well be roundoff errors and not truncation errors
h = surf(log(1e-10 * ones(size(powspec))));
hold on;
g = surf(log(powspec));
hold off;
shading interp; lighting phong;
set(h,'FaceColor',[1 0 0],'FaceAlpha',0.7,'EdgeAlpha', 0);
zlabel('log')
ylabel('|k|');
xlabel('nt');
title('power spectrum estimation');
shg;
pause;

% a quick check of parsevals equation (indicates whether power spectrum
% computation makes sense)
% parseval = zeros(1, Nt);
% for i = 1:Nt
% parseval(i) = (abs(sqrt(sum(powspec(:,i))) - norm(phi(:, i))) );
% end
% plot(parseval);
% title(['parseval, max error = ' num2str(max(parseval))]);
% xlabel('t');
% ylabel('||phi(k)|| - ||phi(x)||');
% shg;
% pause;

% a little movie of the 2d wavefunction
N = sqrt(length(phi(:,1)));
if mod(N,1) == 0 && dim == 2
    for i=1:Nt
        subplot(1,2,1)
        surf(reshape(phi(:,i),N,N))
        title('phi')
        shading interp
        lighting phong
        subplot(1,2,2)
        surf(reshape(psi(:,i),N,N))
        title('psi')
        shading interp
        lighting phong
        shg;
        pause(0.2)
    end
end

%1d movie if data is small enough
if dim == 1 
    for i=1:Nt
        plot(phi(:,i))
        shg;
        pause()
    end
end