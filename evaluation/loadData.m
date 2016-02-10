name = 'compare';
comp = false;

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
if (comp)
    scaling = mass / 0.01;
else
    scaling = 1;
end
tscal = t * scaling; % adjust to compare to karstens time
% t = t + tk(pos); %shifting t
phiAvg = mean(phi);
dphi = dphi / scaling;
dphiAvg = mean(dphi);
psiAvg = mean(psi);
dpsi = dpsi / scaling;
dpsiAvg = mean(dpsi);
rhoAvg = mean(rho);
H = sqrt(rhoAvg / 3);
Nt = length(t);
N = length(phi(:,1));

rhorms = sqrt( (mean(rho.^2) - rhoAvg.^2) ./ rhoAvg.^2 );
rhormsn = sqrt( mean((rho - repmat(rhoAvg,N,1)).^2) ./ rhoAvg.^2 );

rhoalg = (dphi.^2 + (1e-2)^2 * phi.^2) / 2;
rhoAvgalg = mean(rhoalg);
rhormsalg = sqrt( (mean(rhoalg.^2) - rhoAvgalg.^2) ./ rhoAvgalg.^2 );

if dim == 1
    x = linspace(-pi,pi,N+1)';
    x = x(1:end-1);
    
    phik = repmat(phi0k',N,1) + repmat(phi1k',N,1) .* repmat(cos(x),1,length(ak));
    dphik = repmat(dphi0k',N,1) + repmat(dphi1k',N,1) .* repmat(cos(x),1,length(ak));
    
    rhoalgk = (dphik.^2 + (1e-2)^2 * phik.^2) / 2;
    rhoAvgalgk = mean(rhoalgk);
    rhormsalgk = sqrt( (mean(rhoalgk.^2) - rhoAvgalgk.^2) ./ rhoAvgalgk.^2 );
    
    dpsi1 = zeros(N, Nt-2);
    for i = 1 : length(t)-2
        dpsi1(:,i) = (psi(:,i+2) - psi(:,i)) / (t(i+2) - t(i));
    end
    
    dpsi2 = zeros(N, Nt-2);
    dphi2 = zeros(N, Nt-2);
    for i = 1 : length(t)-2
        d1 = t(i+1) - t(i);
        d2 = t(i+2) - t(i+1);
        dpsi2(:,i) = -d2 / (d1 * (d1 + d2)) * psi(:,i) + ...
            (d2 - d1) / (d2 * d1) * psi(:,i+1) + ...
            d1 / (d2 * (d1 + d2)) * psi(:,i+2);
        dphi2(:,i) = -d2 / (d1 * (d1 + d2)) * phi(:,i) + ...
            (d2 - d1) / (d2 * d1) * phi(:,i+1) + ...
            d1 / (d2 * (d1 + d2)) * phi(:,i+2);
    end
    
    dpsialg = 0.5 * dphi .* (phi - repmat(phiAvg,N,1)) - psi .* repmat(H,N,1);
    
    k = [0:N/2-1 0 -N/2+1:-1]';
    xphi = zeros(N, Nt);
    for i = 1:length(t)
        xphi(:,i) = ifft(1i*k.*fft(phi(:,i)));
    end
    
    tmp = fft(xphi.*dphi)./(1i*repmat(k,1,Nt));
    tmp(abs(repmat(k,1,Nt)) < 1e-10) = 0;
    dpsin = ifft(0.5 * tmp - fft(psi) .* repmat(H,N,1));
    
    phifft = fft(phi);
    phi1 = -(phifft(2,:) + phifft(end,:)) / N;

    dphifft = fft(dphi);
    dphi1 = -(dphifft(2,:) + dphifft(end,:)) / N;
    
%     phi0ksp = spline(ak,phi0k,a);
%     phi1ksp = spline(ak,phi1k,a);
%     aksp = spline(ak,ak,a);
end