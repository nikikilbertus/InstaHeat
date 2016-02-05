close all
dim = 1;
name = 'compare';
mass = 0.1;
scaling = 10^3 * mass;

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
t = h5read(name, '/time');
a = h5read(name, '/a');
rho = h5read(name, '/rho');
phi = h5read(name, '/phi');
psi = h5read(name, '/psi');
powspec = h5read(name, '/power_spectrum');
% compute some further properties
%t = t * 6 * 10^3; % adjust to compare to karstens time
t = t * scaling; % adjust to compare to karstens time
phiAvg = mean(phi);
psiAvg = mean(psi);
rhoAvg = mean(rho);
H = sqrt(rhoAvg / 3);
Nt = length(t);
N = length(phi(:,1));

x = linspace(-pi,pi,N+1)';
x = x(1:end-1);

phifft = fft(phi);
phi1 = -(phifft(2,:) + phifft(end,:)) / N;

phi0ksp = spline(tk,phi0k,t);
phi1ksp = spline(tk,phi1k,t);
aksp = spline(tk,ak,t);