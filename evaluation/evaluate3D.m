name = 'compare';

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
dim = h5read(name, '/dimension');
t = h5read(name, '/time');
a = h5read(name, '/a');

rhomean = h5read(name, '/rho_mean');
rhovar = h5read(name, '/rho_variance');

phimean = h5read(name, '/phi_mean');
phivar = h5read(name, '/phi_variance');

dphimean = h5read(name, '/dphi_mean');
dphivar = h5read(name, '/dphi_variance');

psivar = h5read(name, '/psi_variance');
dpsivar = h5read(name, '/dpsi_variance');

powspec = h5read(name, '/power_spectrum');

N = h5read(name, '/gridpoints_internal');
Nout = h5read(name, '/gridpoints_output');

mass = h5read(name, '/mass');

H = sqrt(rhomean / 3);