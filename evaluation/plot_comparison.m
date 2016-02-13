lt = max(t);
st = 800;
s = 1;

disp('a:');
disp(ak(pos) - a(1));
disp('phi0:');
disp(phi0k(pos) - phiAvg(1));
disp('phi1:');
disp(phi1k(pos) - phi1(1));
disp('dphi0:');
disp(dphi0k(pos) - dphiAvg(1));
disp('dphi1:');
disp(dphi1k(pos) - dphi1(1));
disp('psi1:');
disp(psi1k(pos) - psi1(1));
disp('dpsi1:');
disp(dpsi1k(pos) - dpsi1(1));
disp('rhorms:');
disp(rhormsk(pos) - rhorms(1));
disp('k/H:');
disp(koveraH(pos)*ak(pos) - 1/H(1));

%--------------------------by eye comparison phi---------------------------
subplot(4,1,1)
plot(ak(s:posmax), phi0k(s:posmax).*ak(s:posmax).^(3/2), a, phiAvg'.*a.^(3/2));
ylabel('phi0 a^{3/2}')
xlabel('a')

subplot(4,1,2)
plot(ak(s:posmax), dphi0k(s:posmax).*ak(s:posmax).^(3/2), a, dphiAvg'.*a.^(3/2));
ylabel('dphi0 a^{3/2}')
xlabel('a')

subplot(4,1,3)
plot(ak(s:posmax), phi1k(s:posmax).*ak(s:posmax).^(3/2), a, phi1'.*a.^(3/2));
ylabel('phi1 a^{3/2}')
xlabel('a')

subplot(4,1,4)
plot(ak(s:posmax), dphi1k(s:posmax).*ak(s:posmax).^(3/2), a, dphi1'.*a.^(3/2));
ylabel('dphi1 a^{3/2}')
xlabel('a')

shg
pause
figure


%--------------------------errors phi--------------------------------------
subplot(4,1,1)
plot(a, (phiAvg' - phi0ksp).*a.^(3/2));
ylabel('diff: phi0 a^{3/2}')
xlabel('a')

subplot(4,1,2)
plot(a, (dphiAvg' - dphi0ksp).*a.^(3/2));
ylabel('diff: dphi0 a^{3/2}')
xlabel('a')

subplot(4,1,3)
plot(a, (phi1' - phi1ksp).*a.^(3/2));
ylabel('diff: phi1 a^{3/2}')
xlabel('a')

subplot(4,1,4)
plot(a, (dphi1' - dphi1ksp).*a.^(3/2));
ylabel('diff: dphi1 a^{3/2}')
xlabel('a')

shg
pause
figure

%---------------------by eye comparison k/h, rho, psi----------------------
subplot(4,1,1)
plot(ak(s:posmax), koveraH(s:posmax).*ak(s:posmax), a, 1./H);
ylabel('k/H')
xlabel('a')

subplot(4,1,2)
plot(ak(pos:posmax), rhormsk(pos:posmax), a, rhorms)
ylabel('rhorms')
xlabel('a')

subplot(4,1,3)
plot(ak(s:posmax), psi1k(s:posmax).*ak(s:posmax).^(3/2), a, psi1'.*a.^(3/2));
ylabel('psi1 a^{3/2}')
xlabel('a')

subplot(4,1,4)
plot(ak(s:posmax), dpsi1k(s:posmax).*ak(s:posmax).^(3/2), a, dpsi1'.*a.^(3/2));
ylabel('dpsi1 a^{3/2}')
xlabel('a')

shg
pause
figure

%---------------------errors for k/h, rho, psi-----------------------------
subplot(4,1,1)
plot(a, koveraHsp.*a - 1./H');
ylabel('diff: k/H')
xlabel('a')

subplot(4,1,2)
plot(a, rhorms' - rhormsksp)
ylabel('diff: rhorms')
xlabel('a')

subplot(4,1,3)
plot(a, (psi1' - psi1ksp).*a.^(3/2));
ylabel('diff: psi1 a^{3/2}')
xlabel('a')

subplot(4,1,4)
plot(a, (dpsi1' - dpsi1ksp).*a.^(3/2));
ylabel('diff: dpsi1 a^{3/2}')
xlabel('a')

shg
pause

return

figure
subplot(1,2,1)
I = (tk<lt);
II = (t<lt);
plot(t(II),phiAvg(II)'.*a(II).^(3/2),tk(I),phi0k(I).*ak(I).^(3/2))
legend('phi_0 a^{3/2}','phi_0 a^{3/2} kar')
subplot(1,2,2)
I = (tk<st);
II = (t<st);
plot(t(II),phiAvg(II)'.*a(II).^(3/2),tk(I),phi0k(I).*ak(I).^(3/2))
legend('phi_0 a^{3/2}','phi_0 a^{3/2} kar')

figure
subplot(1,2,1)
I = (tk<lt);
II = (t<lt);
plot(t(II),phi1(II)'.*a(II).^(3/2),tk(I),phi1k(I).*ak(I).^(3/2))%,tk(I),dphi1k(I).*ak(I).^(3/2))
legend('phi_1 a^{3/2}','phi_1 a^{3/2} kar')%,'dphi_0 a^{3/2} kar')
subplot(1,2,2)
I = (tk<st);
II = (t<st);
plot(t(II),phi1(II)'.*a(II).^(3/2),tk(I),phi1k(I).*ak(I).^(3/2))
legend('phi_1 a^{3/2}','phi_1 a^{3/2} kar')

figure
subplot(1,2,1)
surf(x,t(t<st),psi(:,t<st)')
title('psi')
xlabel('x')
ylabel('t')
shading interp; lighting phong
subplot(1,2,2)
surf(x,t(t<lt),psi(:,t<lt)')
title('psi')
xlabel('x')
ylabel('t')
shading interp; lighting phong

figure
subplot(1,2,1)
surf(x,t(t<st),rho(:,t<st)')
title('rho')
xlabel('x')
ylabel('t')
shading interp; lighting phong
subplot(1,2,2)
surf(x,t(t<lt),rho(:,t<lt)')
title('rho')
xlabel('x')
ylabel('t')
shading interp; lighting phong

figure
% surf(x,t(t<lt),((rho(:,t<lt) - repmat(rhoAvg(t<lt),N,1)) ./ repmat(max(rho(:,t<lt)),N,1))')
% title('(rho - <rho>) / max(rho)')
surf(x,t(t<lt),((rho(:,t<lt) - repmat(rhoAvg(t<lt),N,1)))')
title('(rho - <rho>)')
xlabel('x')
ylabel('t')
shading interp; lighting phong
