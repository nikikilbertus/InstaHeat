%% evaluate fields by simply cycling through some plots for a first impression

%% setup (user input)
name = 'resolution_short_hardcut/96_16_5e-3';
% one could scale all fields and/or the standard deviations by some function
readDsets;
scal = ones(size(a));
% scal = a.^(3/2);
scalstd = ones(size(a));

%% run
loadDsets;
require(dsetsSummary);
readDsets;
loglog(a,phimax.*scal,a,abs(phimin).*scal,a,abs(phimean).*scal); xlabel('a'); ylabel('\phi'); shg; pause;
loglog(a,phistd.*scalstd); xlabel('a'); ylabel('std \phi'); shg; pause;
loglog(a,dphimax.*scal,a,abs(dphimin).*scal,a,abs(dphimean).*scal); xlabel('a'); ylabel('d\phi'); shg; pause;
loglog(a,dphistd.*scalstd); xlabel('a'); ylabel('std d\phi'); shg; pause;
loglog(a,psimax.*scal,a,abs(psimin).*scal,a,abs(psimean).*scal); xlabel('a'); ylabel('\psi'); shg; pause;
loglog(a,psistd.*scalstd); xlabel('a'); ylabel('std \psi'); shg; pause;
loglog(a,dpsimax.*scal,a,abs(dpsimin).*scal,a,abs(dpsimean).*scal); xlabel('a'); ylabel('d\psi'); shg; pause;
loglog(a,dpsistd.*scalstd); xlabel('a'); ylabel('std d\psi'); shg; pause;
loglog(a,rhomax.*scal,a,abs(rhomin).*scal,a,abs(rhomean).*scal); xlabel('a'); ylabel('\rho'); shg; pause;
loglog(a,rhostd.*scalstd); xlabel('a'); ylabel('std \rho'); shg; pause;
loglog(a,rhorms.*scalstd); xlabel('a'); ylabel('std \rho / |<\rho>|'); shg; pause;
loglog(a,abs(pressuremean).*scal,a,abs(pressuremin).*scal,a,abs(pressuremean).*scal); xlabel('a'); ylabel('p'); shg; pause;
loglog(a,pressurestd.*scalstd); xlabel('a'); ylabel('std p'); shg; pause;