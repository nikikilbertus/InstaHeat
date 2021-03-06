%% DEPRECATED!
plot(t, phiAvg'.*a.^(3/2), tk, phi0k.*ak.^(3/2))
xlabel('t');
ylabel('a^{3/2} * phi_0');
shg
pause
if dim ==1
    plot(t, phi1'.*a.^(3/2), tk, phi1k.*ak.^(3/2))
    xlabel('t');
    ylabel('a^{3/2} * phi_1');
    shg
    pause
end

% interpolate to t
% phi0ksp = spline(tk,phi0k,t);
% phi1ksp = spline(tk,phi1k,t);
% aksp = spline(tk,ak,t);
% subplot(1,2,1)
% plot(t, phiAvg'.*a.^(3/2) - phi0ksp.*aksp.^(3/2))
% xlabel('t');
% ylabel('error in: a^{3/2} * phi_0')
% subplot(1,2,2)
% plot(t, phi1'.*a.^(3/2) - phi1ksp.*aksp.^(3/2))
% xlabel('t');
% ylabel('error in: a^{3/2} * phi_1');
% shg
% pause


% interpolate to tk
% phiAvgsp = spline(t,phiAvg,tk);
% asp = spline(t,a,tk);
% plot(tk, phiAvgsp.*asp.^(3/2) - phi0k.*ak.^(3/2))
% xlabel('t');
% ylabel('error in: a^{3/2} * phi_0');
% shg
% pause

% start with some simple plots (rho over a with a^-4 for reference)
% loglog(a, rhoAvg);
% 'logfit' works only with logfit package at
% http://www.mathworks.com/matlabcentral/fileexchange/29545
% (just put the logfit.m file in the same folder as this one)
% logfit(a, rhoAvg, 'loglog')

% loglog(a, rhoAvg);
% xlabel('a')
% ylabel('rho');
% shg;
% pause;

% H over a with a^-2 for reference
% loglog(a, H);
% logfit(a, H, 'loglog')

% loglog(a, H);
% xlabel('a')
% ylabel('H');
% shg;
% pause;

% average field value and perturbation (with max and min) over time
% subplot(1,3,1)
% plot(t, phiAvg, t, max(phi), t, min(phi));
% title('phi')
% xlabel('t')
% legend('<phi>', 'max', 'min')
% subplot(1,3,2)
% plot(t, rhoAvg, t, max(rho), t, min(rho));
% title('rho')
% xlabel('t')
% legend('<rho>', 'max', 'min')
% subplot(1,3,3)
% plot(t, psiAvg, t, max(psi), t, min(psi));
% title('psi')
% xlabel('t')
% legend('<psi>', 'max', 'min')
% shg;
% pause;
% close all;

% the power spectrum with a reference plane at 10^-10 (everything below
% that might as well be roundoff errors and not truncation errors
% h = surf(log(1e-10 * ones(size(powspec))));
% hold on;
% g = surf(log(powspec));
% hold off;
% shading interp; lighting phong;
% set(h,'FaceColor',[1 0 0],'FaceAlpha',0.7,'EdgeAlpha', 0);
% zlabel('log')
% ylabel('|k|');
% xlabel('nt');
% title('power spectrum estimation');
% shg;
% pause;

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
if dim == 2
    N = sqrt(length(phi(:,1)));
    if mod(N,1) == 0 && dim == 2
        for i=1:Nt
            phiplot = reshape(phi(:,i),N,N);
            psiplot = reshape(psi(:,i),N,N);
            rhoplot = reshape(rho(:,i),N,N);

            subplot(2,3,1)
            surf(phiplot)
            title(['phi, t=' num2str(t(i))])
            shading interp
            lighting phong
            subplot(2,3,4)
            contourf(phiplot)

            subplot(2,3,2)
            surf(rhoplot)
            title('rho')
            shading interp
            lighting phong
            subplot(2,3,5)
            contourf(rhoplot)

            subplot(2,3,3)
            surf(psiplot)
            title('psi')
            shading interp
            lighting phong
            subplot(2,3,6)
            contourf(psiplot)
            shg;
            pause(0.1)
        end
    end
end

%1d movie if data is small enough
if dim == 1
    close all
    Nk = length(phi(:,1));
    x = linspace(-pi,pi,Nk+1)';
    x = x(1:end-1);
    if dim == 1 
        for i=1:Nt
            subplot(1,3,1)
            plot(x,phi(:,i))
            title('phi')
            subplot(1,3,2)
            plot(x,rho(:,i))
            title('rho')
            subplot(1,3,3)
            plot(x,psi(:,i))
            title('psi')
            shg;
            pause(0.1);
        end
    end
end