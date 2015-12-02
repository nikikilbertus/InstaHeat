Nx = 32;
Ny = 32;
Nz = 32;
Ntot = Nx * Ny * Nz;
prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';

name = [prefix 't.txt'];
t = csvread(name);

name = [prefix 'frw_a.txt'];
a = csvread(name);

name = [prefix 'rho.txt'];
rho = csvread(name);

name = [prefix 'pow_spec_000.txt'];
rawPowspec = csvread(name);

name = [prefix 'field_000.txt'];
rawField = csvread(name);

Nt = length(t);
powspec = reshape(rawPowspec, length(rawPowspec)/(Nt-1), Nt-1);
phi = reshape(rawField, Ntot, Nt);
phiAvg = mean(phi);