%% generate and write data to csv file for pgfplots
% modify setup, run init and then choose on of the output sections

%% setup (user input)
% choose (roughly) 'nout' many logarithmically spaced time slices for output
nout = 1500;
name = 'filename';

%% init: get necessary indices (don't change)
readDsets
nouttmp = 0; i = 0;
while (nouttmp < nout) 
    readI = unique(round(logspace(0.1,log10(length(a)),nout+i)));
    nouttmp = length(readI);
    i = i+1;
end

%% output for instability band
require('rhoS','phiS','psiS','mass','spatial_bounds_x','N')
readDsets
L = spatial_bounds_x(2)-spatial_bounds_x(1);
kmingrid = 2*pi/L;
kmaxgrid = sqrt(3)*kmingrid*N(1)/2;
lC = 1./(3*H*mass);
lH = 1./H;
kmin = a/kmingrid;
kmax = a/kmaxgrid;
T = table(a,lH,lC,kmin,kmax);
writetable(T,name);

%% output for evolution overview
require('rhoS','psiS')
readDsets
rhosca = rhomean.*a.^3;
delmaxrho = max(-rhomin,rhomax)./rhomean;
T = table(a,rhorms,rhosca,delmaxrho,stdpsi,maxpsi);
writetable(T,name);

% preview of what this could look like
loglog(a,psistd,a,max(-psimin,psimax),a,rhorms,a,rhomean.*a.^(3),a,max(-rhomin,rhomax)./rhomean,a(((a>4) & (a<150))),0.004*a(((a>4) & (a<150))),'--')