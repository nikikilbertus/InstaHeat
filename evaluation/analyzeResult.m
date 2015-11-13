function [err, energyloss, Tf] = analyzeResult(num, Nx, x, sol, tf, dtf, plotopt)

prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';
name = [prefix 'field_00' int2str(num) '.txt'];
rawField = csvread(name);

Nt = length(rawField) / Nx / 2;
dt = dtf(num);
if ceil(tf / dt) ~= Nt
    error('somethings wrong, check the parameters!')
end
t = linspace(0, (Nt-1) * dt, Nt);

shapedField = reshape(rawField, 2*Nx, Nt);
phi = shapedField(1:Nx, : );

name = [prefix 'a_00' int2str(num) '.txt'];
frwa = csvread(name);

name = [prefix 'energy_00' int2str(num) '.txt'];
rho = csvread(name);
energyloss = (max(rho) - min(rho));

if plotopt == 1
surf(x, t, phi');
title(['numerics, dt = ' num2str(dt)]);
shading interp; lighting phong;
shg; pause();

% surf(x, t, exact);
% title('exact');
% shading interp; lighting phong;
% shg; pause();
% surf(x, t, abs(exact - phi'));
% title('abs difference');
% shading interp; lighting phong;
% shg; pause();

% plot(x, exact(end, : ), x, phi( : , end));
% title(['last timeslice, dt = ' num2str(dt)]);
% shg; pause();

% 
% plot(t,frwa)
% shg
% pause()
% slope = logfit(t, frwa, 'loglog');
% title(['scaling a, dt = ' num2str(dt) ', power = ' num2str(slope)]);
% shg; pause();
% 
% plot(t,rho);
% title('rho');
% xlabel('t');
% ylabel('rho');
% shg;
% pause();

% slope = logfit(t, rho, 'loglog');
% title(['rho, dt = ' num2str(dt) ', diff = ' num2str(energyloss) ', power = ' num2str(slope)]);
% shg; pause();

slope = logfit(frwa, rho, 'loglog');
title(['dt = ' num2str(dt) ', slope = ' num2str(slope)]);
xlabel('a')
ylabel('rho');
shg;
pause()

end

% [X, T] = meshgrid(x, t);
% exact = sol(X, T);
% err = norm(exact(end, : ) - phi( : , end)');

Tf = phi( : , round(floor(tf) / dt));
err = 0;

end