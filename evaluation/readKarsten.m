close all
dim = 1;
name = 'compare';
scaleamp = 1/2;

name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten_compare/' name '.dat'];

raw = importdata(name);

tk = raw(:,1);
ak = raw(:,2);
phi0k = raw(:,3) .* ak.^(-3/2);
dphi0k = raw(:,4) .* ak.^(-3/2);
phi1k = raw(:,5) .* ak.^(-3/2) * scaleamp;
dphi1k = raw(:,6) .* ak.^(-3/2) * scaleamp;
delrhok = raw(:,8);

plot(tk, phi0k.*ak.^(3/2))
xlabel('t');
ylabel('a^{3/2} * phi_0');
shg;
pause;

plot(tk, phi1k.*ak.^(3/2))
xlabel('t');
ylabel('a^{3/2} * phi_1');
shg;
pause;

% plot(tk, delrhok)
% xlabel('t');
% ylabel('delta rho');
% shg;
% pause;

clear raw