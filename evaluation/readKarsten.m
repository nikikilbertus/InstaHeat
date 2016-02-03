close all
dim = 1;
Nk = 512;
name = 'compare';

name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten_compare/' name '.dat'];

raw = importdata(name);

tk = raw(:,1);
ak = raw(:,2);
phi0k = raw(:,3) .* ak.^(-3/2);
dphi0k = raw(:,4) .* ak.^(-3/2);
phi1k = raw(:,5) .* ak.^(-3/2);
dphi1k = raw(:,6) .* ak.^(-3/2);
delrhok = raw(:,8);

x = linspace(-pi,pi,Nk+1)';
x = x(1:end-1);

phik = zeros(Nk,length(tk));
dphik = zeros(Nk,length(tk));
for i = 1:length(tk)
    phik(:,i) = phi0k(i) + phi1k(i) * cos(x);
    dphik(:,i) = dphi0k(i) + dphi1k(i) * cos(x);
end

plot(tk, phi0k.*ak.^(3/2))
xlabel('t');
ylabel('a^(3/2) * phi0');
shg;
pause;

plot(tk, dphi0k.*ak.^(3/2))
xlabel('t');
ylabel('a^(3/2) * phi0');
shg;
pause;

clear raw