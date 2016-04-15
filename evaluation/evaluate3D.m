% name = 'comp_karsten_5_simp2';
interp = false;
newoutput = true;
cstr = false;

% loading the data, replace 'name' with the path where you stored the .h5
% file from the simulation
name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name '.h5'];

% built in functions for reading makes it easy
dim = h5read(name, '/dimension');
t = h5read(name, '/time');
a = h5read(name, '/a');
N = h5read(name, '/gridpoints_internal');
Nout = h5read(name, '/gridpoints_output');
mass = h5read(name, '/mass');
strides = h5read(name,'/strides_space');
timestride = h5read(name,'/strides_time');
tols = h5read(name, '/tolerances');

if (newoutput)
    powspec  = h5read(name, '/phi_power_spectrum');
    phismry  = h5read(name, '/phi_summary');
    dphismry = h5read(name, '/dphi_summary');
    psismry  = h5read(name, '/psi_summary');
    dpsismry = h5read(name, '/dpsi_summary');
    rhosmry  = h5read(name, '/rho_summary');
    phimean  = phismry(1,:);
    phivar   = phismry(2,:);
    phimin   = phismry(3,:);
    phimax   = phismry(4,:);
    dphimean = dphismry(1,:);
    dphivar  = dphismry(2,:);
    dphimin  = dphismry(3,:);
    dphimax  = dphismry(4,:);
    psimean  = psismry(1,:);
    psivar   = psismry(2,:);
    psimin   = psismry(3,:);
    psimax   = psismry(4,:);
    dpsimean = dpsismry(1,:);
    dpsivar  = dpsismry(2,:);
    dpsimin  = dpsismry(3,:);
    dpsimax  = dpsismry(4,:);
    rhomean  = rhosmry(1,:);
    rhovar   = rhosmry(2,:);
    rhomin   = rhosmry(3,:);
    rhomax   = rhosmry(4,:);
    if (cstr)
       cstr = h5read(name, '/constraints'); 
    end
else
    powspec = h5read(name, '/power_spectrum');
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
end
H = sqrt(rhomean / 3);
rhorms = sqrt(rhovar ./ rhomean.^2);

if interp
[phipks, phipkpos] = findpeaks(phimean);
phienv = spline(a(phipkpos), phipks, a);
end
if interp
[phipks, phipkpos] = findpeaks(phivar);
phivarenv = spline(a(phipkpos), phipks, a);
phirms = sqrt(phivarenv ./ phienv.^2);
end

if interp
[dphipks, dphipkpos] = findpeaks(dphimean);
dphienv = spline(a(dphipkpos), dphipks, a);
end
if interp
[dphipks, dphipkpos] = findpeaks(dphivar);
dphivarenv = spline(a(dphipkpos), dphipks, a);
dphirms = sqrt(dphivarenv ./ dphienv.^2);
end