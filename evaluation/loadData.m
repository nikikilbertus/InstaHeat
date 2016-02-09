name = 'compare_2';

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
dim = h5read(name, '/dimension');
t = h5read(name, '/time');
a = h5read(name, '/a');
rho = h5read(name, '/rho');
phi = h5read(name, '/phi');
psi = h5read(name, '/psi');
powspec = h5read(name, '/power_spectrum');
mass = h5read(name, '/mass');
% compute some further properties
scaling = mass / 0.01;
t = t * scaling; % adjust to compare to karstens time
% t = t + tk(pos); %shifting t
phiAvg = mean(phi);
psiAvg = mean(psi);
rhoAvg = mean(rho);
H = sqrt(rhoAvg / 3);
Nt = length(t);
N = length(phi(:,1));

if dim == 1
    x = linspace(-pi,pi,N+1)';
    x = x(1:end-1);

    phifft = fft(phi);
    phi1 = -(phifft(2,:) + phifft(end,:)) / N;

%     phi0ksp = spline(ak,phi0k,a);
%     phi1ksp = spline(ak,phi1k,a);
%     aksp = spline(ak,ak,a);
end