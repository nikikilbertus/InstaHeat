%% evaluate gravitational power spectrum from single run

%% setup (user input)
% the filename
name = 'masses96/96_0.005';
% the number of timeslices for which to show the power spectrum
nslices = 60;
fieldname = 'psi';

%% run
% initialization
require([fieldname 'ps'],'rhoS','N','spatial_bounds_x');
readDsets;
field = psips;
nbins = size(field,2);
N = N(1);
L = spatial_bounds_x(2)-spatial_bounds_x(1);
kmin = 2*pi/L; kmax = sqrt(3)*kmin*N/2;
k = (1:nbins)*kmax/nbins;
% check for negative and nan values (just making sure)
I = (field<0);
if any(I(:))
    display('found negative entries, set to 0');
    field(I) = 0;
end
I = isnan(field);
if any(I(:))
    display('found nan entries, set to 0');
    field(I) = 0;
end
% body
skip = int64(floor(length(a)/nslices));
J = jet; % colormap
subplot(1,2,2);
loglog(a,rhorms); hold on;
% range = 1:skip:length(a);
I = getLogIndices(length(a),nslices);
range = 1:length(a);
range = range(I);
for i = [range; 1:length(range)]
    subplot(1,2,1)
    loglog(k,field(i(1),:),'color',J(i(2),:));
    xlabel('|k|'); ylabel(['d \Omega_{' fieldname '} / d ln(k)']); hold on;
    subplot(1,2,2)
    loglog(a(i(1)),rhorms(i(1)),'x','color',J(i(2),:));
    xlabel('a'); ylabel('std(\rho) / |<\rho>|');
end
% axis([-inf inf 1e-20 inf]);
hold off; shg;