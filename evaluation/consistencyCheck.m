%% DEPRECATED!
tmpres1 = xxpsi ./ repmat(a.^2',N,1) - 3 * repmat(H.^2,N,1) .* psi - ...
    3 * repmat(H,N,1) .* dpsi - 0.5 * (rho - repmat(rhoAvg,N,1));

figure
surf(t,x,tmpres1)
% plot(t,mean(tmpres1))
xlabel('t')
ylabel('x')
title('00 equation with dpsi from code')
shading interp
shg
pause

tmpres2 = xxpsi ./ repmat(a.^2',N,1) - 3 * repmat(H.^2,N,1) .* psi - ...
    3 * repmat(H,N,1) .* dpsipad - 0.5 * (rho - repmat(rhoAvg,N,1));

figure
surf(t(2:end-1),x,tmpres2(:,2:end-1))
% plot(t(2:end-1),mean(tmpres2(:,2:end-1)))
xlabel('t')
ylabel('x')
title('00 equation with dpsi from findiff')
shading interp
shg
pause

figure
surf(t,x,ifft(repmat(1i*k,1,Nt).*fft(dpsi + repmat(H,N,1) .* psi)) - 0.5 * xphi .* dphi)
% plot(t,mean(ifft(repmat(1i*k,1,Nt).*fft(dpsi + repmat(H,N,1) .* psi)) - 0.5 * xphi .* dphi))
xlabel('t')
ylabel('x')
title('0i equation with dpsi from code')
shading interp
shg
pause

figure
tmp = ifft(repmat(1i*k,1,Nt).*fft(dpsipad + repmat(H,N,1) .* psi)) - 0.5 * xphi .* dphi;
surf(t(2:end-1),x,tmp(:,2:end-1))
% plot(t(2:end-1),mean(tmp(:,2:end-1)))
xlabel('t')
ylabel('x')
title('0i equation with dpsi from findiff')
shading interp
shg
pause

figure
% surf(t,x,econ)
plot(t,mean(econ))
xlabel('t')
% ylabel('x')
title('nabla_a T^{a0} with ddphi from eom')
shading interp
shg
pause

figure
% surf(t(2:end-1),x,econfd(:,2:end-1))
plot(t(2:end-1),mean(econfd(:,2:end-1)))
xlabel('t')
% ylabel('x')
title('nabla_a T^{a0} with ddphi & dpsi from findiff')
shading interp
shg


% clear
% readKarsten
% close all
% econs = {}; econsfd={}; ts = {};
% for kk = 4:2:18
% name = ['compare_' num2str(kk)]
% loadData
% econs{int16(kk/2-1)} = econ;
% econsfd{int16(kk/2-1)} = econfd;
% ts{int16(kk/2-1)} = t;
% end

% figure
% hold on
% for i = 8:-1:1
% plot(ts{i},mean(econs{i}))
% end
% hold off
% shg