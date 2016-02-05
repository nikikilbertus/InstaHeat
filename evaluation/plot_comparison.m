lt = max(t);
st = 800;

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
