% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/run.h5';

% built in functions for reading makes it easy
t = h5read(name, '/time');
a = h5read(name, '/a');
rho = h5read(name, '/rho');
phi = h5read(name, '/phi');
powspec = h5read(name, '/power_spectrum');
% compute some further properties
H = sqrt(rho / 3);
phiAvg = mean(phi);
Nt = length(t);

% start with some simple plots (rho over a with a^-4 for reference)
loglog(a, rho, a, a.^(-4) * rho(1));
% 'logfit' works only with logfit package at
% http://www.mathworks.com/matlabcentral/fileexchange/29545
% (just put the logfit.m file in the same folder as this one)
% logfit(a, rho, 'loglog')
xlabel('a')
ylabel('rho');
shg;
pause;

% H over a with a^-2 for reference
loglog(a, H, a, a.^(-2) * H(1));
% logfit(a, H, 'loglog')
xlabel('a')
ylabel('H');
shg;
pause;

% average field value (with max and min) over time
plot(t, phiAvg, t, max(phi), t, min(phi));
xlabel('t')
legend('<phi>', 'max', 'min')
shg;
pause;

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
shg;
pause;

% a quick check of parsevals equation (indicates whether power spectrum
% computation makes sense)
parseval = zeros(1, Nt);
for i = 1:Nt
parseval(i) = (abs(sqrt(sum(powspec(:,i))) - norm(phi(:, i))) );
end
plot(parseval);
title(['parseval, max error = ' num2str(max(parseval))]);
xlabel('t');
ylabel('||phi(k)|| - ||phi(x)||');
shg;
pause;