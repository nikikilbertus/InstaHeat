name = 'test_hyp_6000_scal-2';
comp = true;
karstenpsi = true;
loadtimes = true;

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
dim = h5read(name, '/dimension');
t = h5read(name, '/time')';
a = h5read(name, '/a')';
rho = h5read(name, '/rho');
phi = h5read(name, '/phi');
dphi = h5read(name, '/dphi');
psi = h5read(name, '/psi');
dpsi = h5read(name, '/dpsi');
powspec = h5read(name, '/power_spectrum');
mass = h5read(name, '/mass');
tols = h5read(name, '/tolerances');

Nt = length(t);
N = length(phi(:,1));

phimean = h5read(name, '/phi_mean')';
phivar = h5read(name, '/phi_variance')';
dphimean = h5read(name, '/dphi_mean')';
dphivar = h5read(name, '/dphi_variance')';
psimean = h5read(name, '/psi_mean')';
psivar = h5read(name, '/psi_variance')';
dpsimean = h5read(name, '/dpsi_mean')';
dpsivar = h5read(name, '/dpsi_variance')';
rhomean = h5read(name, '/rho_mean')';
rhovar = h5read(name, '/rho_variance')';

H = sqrt(rhomean / 3);

rhorms = sqrt( rhovar ) ./ rhomean;
rhormsalt1   = sqrt( (mean(rho.^2) - rhomean.^2) ./ rhomean.^2 );
rhormsalt2 = sqrt( mean(rho.^2) - rhomean.^2 ) ./ rhomean;
rhormsalt3  = sqrt( mean((rho - repmat(rhomean,N,1)).^2) ./ rhomean.^2 );

% compute some further properties
if (comp)
    scaling = mass / massk;
else
    scaling = 1;
end

if max(ak) > max(a)
    [~, pos] = min((ak - min(a)).^2);
    [~, posmax] = min((ak - max(a)).^2);
    posmax = posmax + min(10, length(ak) - pos);
else
    [~, pos] = min((ak - min(a)).^2);
    posmax = length(ak);
end

tunscaled = t;
t = tunscaled * scaling; % adjust to compare to karstens time
% t = t + tk(pos); %shifting t
dphi = dphi / scaling;
dphimean = dphimean / scaling;
dphivar = dphivar / scaling^2;

dpsi = dpsi / scaling;
dpsimean = dpsimean / scaling;
dpsivar = dpsivar / scaling^2;

% rhoalg = (dphi.^2 + (massk)^2 * phi.^2) / 2;
% rhomeanalg = mean(rhoalg);
% rhormsalg = sqrt( (mean(rhoalg.^2) - rhomeanalg.^2) ./ rhomeanalg.^2 );

if dim == 1
    x = linspace(-pi,pi,N+1)';
    x = x(1:end-1);
    
    k = [0:N/2-1 0 -N/2+1:-1]';
    klap = k;
    klap(round(N/2) + 1) = N/2;
    xphi = ifft(repmat(1i*k,1,Nt).*fft(phi));
    xxphi = ifft(repmat(-klap.^2,1,Nt).*fft(phi));
    xpsi = ifft(repmat(1i*k,1,Nt).*fft(psi));
    xxpsi = ifft(repmat(-klap.^2,1,Nt).*fft(psi));
    
    amat = repmat(a,N,1);
    Hmat = repmat(H,N,1);
    V = mass^2 * phi.^2 / 2;
    Vprime = mass^2 * phi;
    ddphi = (1 + 4*psi).*xxphi ./ amat.^2 - (1 + 2*psi).*V + 4*dpsi.*dphi - 3*Hmat.*dphi;
    
%     phik = repmat(phi0k',N,1) + repmat(phi1k',N,1) .* repmat(cos(x),1,length(ak));
%     dphik = repmat(dphi0k',N,1) + repmat(dphi1k',N,1) .* repmat(cos(x),1,length(ak));
    
%     rhoalgk = (dphik.^2 + (massk)^2 * phik.^2) / 2;
%     rhomeanalgk = mean(rhoalgk);
%     rhormsalgk = sqrt( (mean(rhoalgk.^2) - rhomeanalgk.^2) ./ rhomeanalgk.^2 );
    
    dpsifd1 = zeros(N, Nt-2);
    dpsifd2 = zeros(N, Nt-2);
    dphifd2 = zeros(N, Nt-2);
    ddphifd2 = zeros(N, Nt-2);
    for i = 1 : length(t)-2
        dpsifd1(:,i) = (psi(:,i+2) - psi(:,i)) / (t(i+2) - t(i));
        
        d1 = t(i+1) - t(i);
        d2 = t(i+2) - t(i+1);
        dpsifd2(:,i) = -d2 / (d1 * (d1 + d2)) * psi(:,i) + ...
            (d2 - d1) / (d2 * d1) * psi(:,i+1) + ...
            d1 / (d2 * (d1 + d2)) * psi(:,i+2);
        dphifd2(:,i) = -d2 / (d1 * (d1 + d2)) * phi(:,i) + ...
            (d2 - d1) / (d2 * d1) * phi(:,i+1) + ...
            d1 / (d2 * (d1 + d2)) * phi(:,i+2);
                
        del = (t(i+1) - t(i+2)) * t(i)^2 + ...
              (t(i+2) - t(i)) * t(i+1)^2 + ... 
              (t(i) - t(i+1)) * t(i+2)^2;
          
        ddphifd2(:,i) = 2 * ( (t(i+1) - t(i+2)) * phi(:,i) + ...
            (t(i+2) - t(i)) * phi(:,i+1) + ...
            (t(i) - t(i+1)) * phi(:,i+2) ) / del;
    end
    
    dphipad = zeros(N,Nt);
    dphipad(:,2:end-1) = dphifd2;
    
    ddphipad = zeros(N,Nt);
    ddphipad(:,2:end-1) = ddphifd2;
    
    dpsialg = 0.5 * dphi .* (phi - repmat(phimean,N,1)) - psi .* repmat(H,N,1);
    dpsipad = zeros(N,Nt);
    dpsipad(:,2:end-1) = dpsifd2;
    
    tmp = fft(xphi.*dphi)./(1i*repmat(k,1,Nt));
    tmp(abs(repmat(k,1,Nt)) < 1e-18) = 0;
    dpsin = ifft(0.5 * tmp - fft(psi) .* repmat(H,N,1));
    
    phifft = fft(phi);
    phi1 = -(phifft(2,:) + phifft(end,:)) / N;
    dphifft = fft(dphi);
    dphi1 = -(dphifft(2,:) + dphifft(end,:)) / N;
    
    psifft = fft(psi);
    psi1 = -(psifft(2,:) + psifft(end,:)) / N;
    dpsifft = fft(dpsi);
    dpsi1 = -(dpsifft(2,:) + dpsifft(end,:)) / N;
    
    econ = dphi.*(amat.^2.*(4.*dphi.*dpsi - 3.*dphi.*Hmat.*(1 - 4.*psi) + ...
    ddphi.*(-1 + 4.*psi) - (1 - 2.*psi).*Vprime) + xxphi) ./ amat.^2;

    econfd = dphipad.*(amat.^2.*(4.*dphipad.*dpsipad - 3.*dphipad.*Hmat.*(1 - 4.*psi) + ...
    ddphipad.*(-1 + 4.*psi) - (1 - 2.*psi).*Vprime) + xxphi) ./ amat.^2;
    
    if max(ak) > max(a)
        I = (ak > 0.8*ak(pos));
%         I = (ak > 1);
        phi0ksp = spline(ak(I),phi0k(I),a);
        phi1ksp = spline(ak(I),phi1k(I),a);
        dphi0ksp = spline(ak(I),dphi0k(I),a);
        dphi1ksp = spline(ak(I),dphi1k(I),a);
        if (karstenpsi)
            psi1ksp = spline(ak(I),psi1k(I),a);
            dpsi1ksp = spline(ak(I),dpsi1k(I),a);
            koveraHsp = spline(ak(I),koveraH(I),a);
            rhormsksp = spline(ak(I),rhormsk(I),a);
        end
    else
        I = (ak >= ak(pos));
        phi0sp = spline(a,phimean,ak(I));
        phi1sp = spline(a,phi1,ak(I));
        dphi0sp = spline(a,dphimean,ak(I));
        dphi1sp = spline(a,dphi1,ak(I));
        if (karstenpsi)
            psi1sp = spline(a,psi1,ak(I));
            dpsi1sp = spline(a,dpsi1,ak(I));
            koverHsp = spline(a,1./H,ak(I));
            rhormssp = spline(a,rhorms,ak(I));
        end
    end
    
    check = xxpsi ./ amat.^2 - ...
            3 * Hmat.^2 .* psi - ...
            3 * Hmat .* dpsipad - ...
            0.5 * (rho - repmat(rhomean, N, 1));
        
    dpsicheck = xxpsi ./ (3 * amat.^2 .* Hmat) - ...
                Hmat .* psi - ...
                0.5 * (rho - repmat(rhomean, N, 1)) ./ (3 * Hmat);
end

if loadtimes
    runtime = h5read(name, '/runtime_stepper');
    try
    steps = h5read(name, '/steps_total');
    stepsok = h5read(name, '/steps_ok');
    stepsbad = h5read(name, '/steps_bad');
    catch me
        %
    end
end