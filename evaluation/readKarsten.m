dim = 1;
namepre = 'compare_psi';
scaleamp = 1e3;
massk = 1e-2;

name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten_compare/' namepre '.dat'];

raw = importdata(name);

if strcmp(namepre, 'compare')
    tk = raw(:,1);
    ak = raw(:,2);
    phi0k = raw(:,3) .* ak.^(-3/2);
    dphi0k = raw(:,4) .* ak.^(-3/2);
    phi1k = raw(:,5) .* ak.^(-3/2) * scaleamp;
    dphi1k = raw(:,6) .* ak.^(-3/2) * scaleamp;
    delrhok = raw(:,8);
else
    ak = raw(:,1);
    phi0k = raw(:,2);
    dphi0k = raw(:,3);
    phi1k = raw(:,4) * scaleamp;
    dphi1k = raw(:,5) * scaleamp;
    koveraH = raw(:,6);
    rhormsk = raw(:,7) / sqrt(2 * pi) * scaleamp;
    psi1k = raw(:,8) * scaleamp;
    dpsi1k = raw(:,9) * scaleamp;
end

figure
plot(ak, phi0k.*ak.^(3/2))
xlabel('a');
ylabel('a^{3/2} * phi_0');
shg;
pause;

plot(ak, phi1k.*ak.^(3/2))
xlabel('a');
ylabel('a^{3/2} * phi_1');
shg;
pause;

plot(ak, psi1k.*ak.^(3/2))
xlabel('a');
ylabel('a^{3/2} * psi');
shg;
pause;

plot(ak, dpsi1k.*ak.^(3/2))
xlabel('a');
ylabel('a^{3/2} * dpsi');
shg;
pause;

% plot(ak, delrhok)
% xlabel('a');
% ylabel('delta rho');
% shg;
% pause;

clear raw