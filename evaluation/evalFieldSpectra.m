%% evaluate the field power spectra from a single run

%% setup (user input)
% the filename
name = 'res2d_long/128_32_3e-5';
% on how many time slices do we want to show the power spectra in the video
nslices = 30;

%% run
% initialization
loadDsets;
require(dsetsSpectra,'rhoS','N','mass','spatial_bounds_x');
readDsets;
nbins = size(phips,2);
N = N(1);
L = spatial_bounds_x(2)-spatial_bounds_x(1);
kmin = 2*pi/L; kmax = sqrt(3)*kmin*N/2;
k = linspace(kmin,kmax,nbins);
% TODO what mass do I need here?
lc = 1./sqrt(3*H*mass);
% body
for ii = [1 logspace(1,log10(length(a)),nslices)]
    i = int64(floor(ii));
    subplot(2,2,1)
    loglog(k,phips(i,1:nbins)); xlabel('k'); ylabel('power'); hold on
    ktmp = (1/lc(i))*a(i);
    kmax = H(i)*a(i);
    lowb = min(phips(i,1:nbins))+1e-13;
    higb = max(phips(i,1:nbins));
    loglog([ktmp ktmp],[lowb higb]); loglog([kmax kmax],[lowb higb]); hold off;
    title(['\phi, a = ' num2str(a(i)) ' / ' num2str(a(end))]);
    subplot(2,2,2)
    plot(a,rhorms,a(i),rhorms(i),'or'); xlabel('a'); ylabel('std \rho / <\rho>');
    shg;
    subplot(2,2,3)
    loglog(k,rhops(i,1:nbins)); xlabel('k'); ylabel('power'); title('\rho');
    hold on
    ktmp = H(1);
    lowb = min(rhops(i,1:nbins))+1e-10;
    higb = max(rhops(i,1:nbins));
    loglog([ktmp ktmp],[lowb higb]);
    hold off
    subplot(2,2,4)
    if exist('psips','var')
        loglog(k,psips(i,1:nbins)); xlabel('k'); ylabel('power'); title('\psi');
    end
    if i == 1
        Hk = H(i);
        lowb = min(rhops(i,1:nbins))+1e-13;
        higb = max(rhops(i,1:nbins));
        subplot(2,2,3); hold on; plot([Hk Hk], [lowb higb]); hold off;
        if exist('psips','var')
            lowb = min(psips(i,1:nbins))+1e-13;
            higb = max(psips(i,1:nbins));
            subplot(2,2,4); hold on; plot([Hk Hk], [lowb higb]); hold off;
        end
        pause;
    else
        pause(0.05);
    end
end