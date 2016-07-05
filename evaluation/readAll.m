%% read all available datasets of one run

%% setup (user input)
% the filename
name = 'test32';

%% run
loadDsets
dsets = dsetsAll;
readDsets;
Nt = length(a);
phi = reshape(phi,Nt,N(1),N(2),N(3));
dphi = reshape(dphi,Nt,N(1),N(2),N(3));
psi = reshape(psi,Nt,N(1),N(2),N(3));
dpsi = reshape(dpsi,Nt,N(1),N(2),N(3));
rho = reshape(rho,Nt,N(1),N(2),N(3));