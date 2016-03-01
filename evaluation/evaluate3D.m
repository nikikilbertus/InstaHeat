name = 'comp3d_128_5e5_old';
interp = true;

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
dim = h5read(name, '/dimension');
t = h5read(name, '/time');
a = h5read(name, '/a');

rhomean = h5read(name, '/rho_mean');
rhovar = h5read(name, '/rho_variance');
rhorms = sqrt(rhovar ./ rhomean.^2);

phimean = h5read(name, '/phi_mean');
if interp
[phipks, phipkpos] = findpeaks(phimean);
phienv = spline(a(phipkpos), phipks, a);
end
phivar = h5read(name, '/phi_variance');
if interp
[phipks, phipkpos] = findpeaks(phivar);
phivarenv = spline(a(phipkpos), phipks, a);
phirms = sqrt(phivarenv ./ phienv.^2);
end

dphimean = h5read(name, '/dphi_mean');
if interp
[dphipks, dphipkpos] = findpeaks(dphimean);
dphienv = spline(a(dphipkpos), dphipks, a);
end
dphivar = h5read(name, '/dphi_variance');
if interp
[dphipks, dphipkpos] = findpeaks(dphivar);
dphivarenv = spline(a(dphipkpos), dphipks, a);
dphirms = sqrt(dphivar ./ dphienv.^2);
end

psivar = h5read(name, '/psi_variance');
dpsivar = h5read(name, '/dpsi_variance');

powspec = h5read(name, '/power_spectrum');

N = h5read(name, '/gridpoints_internal');
Nout = h5read(name, '/gridpoints_output');

mass = h5read(name, '/mass');

strides = h5read(name,'/strides_space');
timestride = h5read(name,'/strides_time');
tols = h5read(name, '/tolerances');

H = sqrt(rhomean / 3);