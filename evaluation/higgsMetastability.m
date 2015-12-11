prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';
name = [prefix 'run.h5'];

t = h5read(name, '/time');

a = h5read(name, '/a');

rho = h5read(name, '/rho');

H = sqrt(rho / 3);

phi = h5read(name, '/phi');

powspec = h5read(name, '/power_spectrum');

Nt = length(t);

phiAvg = mean(phi);

loglog(a, rho);
% logfit(a, rho, 'loglog')
xlabel('a')
ylabel('rho');
shg;
pause;

loglog(a, H);
% logfit(a, H, 'loglog')
xlabel('a')
ylabel('H');
shg;
pause;

plot(t, phiAvg, t, max(phi), t, min(phi));
xlabel('t')
ylabel('<phi>');
shg;
pause;

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

% surf(log(powspec + 1e-10));
% shading interp; lighting phong;
% zlabel('lin')
% ylabel('|k|');
% xlabel('nt');
% shg;
% pause;

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