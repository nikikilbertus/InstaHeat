name = 'compare';

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
dim = h5read(name, '/dimension');
t = h5read(name, '/time');
a = h5read(name, '/a');
rho = h5read(name, '/rho');
phi = h5read(name, '/phi');
dphi = h5read(name, '/dphi');
psi = h5read(name, '/psi');
dpsi = h5read(name, '/dpsi');
powspec = h5read(name, '/power_spectrum');
mass = h5read(name, '/mass');

% compute some further properties
scaling = mass / 0.01;
%t = t * scaling; % adjust to compare to karstens time
% t = t + tk(pos); %shifting t
phiAvg = mean(phi);
dphiAvg = mean(dphi) / scaling;
psiAvg = mean(psi);
dpsiAvg = mean(dpsi) / scaling;
rhoAvg = mean(rho);
H = sqrt(rhoAvg / 3);
Nt = length(t);
N = length(phi(:,1));

rhorms = sqrt( (mean(rho.^2) - rhoAvg.^2) ./ rhoAvg.^2 );
rhormsn = sqrt( mean((rho - repmat(rhoAvg,N,1)).^2) ./ rhoAvg.^2 );

if dim == 1
    x = linspace(-pi,pi,N+1)';
    x = x(1:end-1);
    
    dpsi1 = zeros(N, Nt-1);
    for i = 1 : length(t)-1
        dpsi1(:,i) = (psi(:,i+1) - psi(:,i)) / (t(i+1) - t(i));
    end
    
    dpsin = 0.5 * dphi .* (phi - repmat(phiAvg,N,1)) - psi .* repmat(H,N,1);
    
    k = [0:N/2-1 0 -N/2+1:-1]';
    xphi = zeros(N, Nt);
    for i = 1:length(t)
        xphi(:,i) = ifft(1i*k.*fft(phi(:,i)));
    end
    
    tmp = fft(xphi.*dphi)./(1i*repmat(k,1,Nt));
    tmp(abs(repmat(k,1,Nt)) < 1e-10) = 0;
    dpsinn = ifft(0.5 * tmp - fft(psi) .* repmat(H,N,1));
    
    phifft = fft(phi);
    phi1 = -(phifft(2,:) + phifft(end,:)) / N;

    dphifft = fft(dphi);
    dphi1 = -(dphifft(2,:) + dphifft(end,:)) / N;
    
%     phi0ksp = spline(ak,phi0k,a);
%     phi1ksp = spline(ak,phi1k,a);
%     aksp = spline(ak,ak,a);
end