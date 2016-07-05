using PyPlot
using HDF5

dim = 2
name = "run"

name = join(["/Users/niki/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/" name ".h5"], "")

t = h5read(name, "/time");
a = h5read(name, "/a");
rho = h5read(name, "/rho");
phi = h5read(name, "/phi");
psi = h5read(name, "/psi");
powspec = h5read(name, "/power_spectrum");

H = sqrt(rho / 3)
phiAvg = mean(phi, 1)
psiAvg = mean(psi, 1)
Nt = length(t)

loglog(a, rho, a, a.^(-4) * rho[1])
xlabel("a")
ylabel(L"\rho")
legend(["data", "reference: a^-4"])
readline(STDIN)
cla()

subplot(1,2,1)
plot(t, squeeze(phiAvg,1), t, squeeze(maximum(phi,1),1), t, squeeze(minimum(phi,1),1))
title(L"\phi")
xlabel("t")
legend([L"<\phi>", "max", "min"])
subplot(1,2,2)
plot(t, squeeze(psiAvg,1), t, squeeze(maximum(psi,1),1), t, squeeze(minimum(psi,1),1))
title(L"\psi")
xlabel("t")
legend([L"<\psi>", "max", "min"])
readline(STDIN)
cla()

surf(log(powspec))
readline(STDIN)
cla()

N = sqrt(length(phi[:,1]))
if isinteger(N) && dim == 2
    N = Int64(N)
    for i = 1:Nt
        phiplot = reshape(phi[:,i], N, N)
        psiplot = reshape(psi[:,i], N, N)
        subplot(2,2,1)
        cla()
        surf(phiplot)
        title(L"\phi")
        subplot(2,2,2)
        cla()
        surf(psiplot)
        title(L"\psi")
        subplot(2,2,3)
        cla()
        contourf(phiplot)
        subplot(2,2,4)
        cla()
        contourf(psiplot)
        pause(0.05)
    end
end
