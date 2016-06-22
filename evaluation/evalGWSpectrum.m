%% evaluate gravitational power spectrum from single run

%% setup
% the filename
name = 'cutoff/96_16_5e-3';
% the number of timeslices for which to show the power spectrum
nslices = 60;

%% run
% initialization
require('gwps','rhoS','N','spatial_bounds_x');
readDsets;
nbins = size(gwps,2);
N = N(1);
L = spatial_bounds_x(2)-spatial_bounds_x(1);
kmin = 2*pi/L; kmax = sqrt(3)*kmin*N/2;
k = (1:nbins)*kmax/nbins;
% check for negative and nan values (just making sure)
I = (gwps<0);
if any(I(:))
    display('found negative entries, set to 0');
    gwps(I) = 0;
end
I = isnan(gwps);
if any(I(:))
    display('found nan entries, set to 0');
    gwps(I) = 0;
end
% body
skip = int64(floor(length(a)/nslices));
J = jet; % colormap
subplot(1,2,2);
loglog(a,rhorms); hold on;
range = 1:skip:length(a);
for i = [range; 1:length(range)]
    subplot(1,2,1)
    loglog(k,gwps(i(1),:),'color',J(i(2),:));
    xlabel('|k|'); ylabel('d \Omega_{gw} / d ln(k)'); hold on;
    subplot(1,2,2)
    loglog(a(i(1)),rhorms(i(1)),'x','color',J(i(2),:));
    xlabel('a'); ylabel('std(\rho) / |<\rho>|');
end
hold off; shg;