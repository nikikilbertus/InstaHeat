%% show the instability band and the modes evolved in a singel run

%% setup (user input)
name = 'res2daggr/256_1e-4_8'

%% output for instability band
require('rhoS','mass','spatial_bounds_x','spatial_bounds_y','spatial_bounds_z','N')
readDsets
Lx = spatial_bounds_x(2)-spatial_bounds_x(1); kx = 2*pi/Lx;
Ly = spatial_bounds_y(2)-spatial_bounds_y(1); ky = 2*pi/Ly;
Lz = spatial_bounds_z(2)-spatial_bounds_z(1); kz = 2*pi/Lz;
kmingrid = min([kx,ky,kz]);
kmaxgrid = sqrt((kx*floor(N(1)/2))^2+(ky*floor(N(2)/2))^2+(kz*floor(N(2)/2)));
lC = sqrt(1./(3*H*mass));
lH = 1./H;
kmin = a/kmingrid;
kmax = a/kmaxgrid;

loglog(a,lH,a,lC,a,kmin,a,kmax);
legend('l_H','l_C','k_{min}','k_{max}'); shg;