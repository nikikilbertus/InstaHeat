%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume visualization video
Nt = length(a);
rho = h5read(name, '/rho');
rho = reshape(rho, N(1), N(2), N(3), Nt);

x = linspace(0,1,N(1)+1)';
x = x(1:end-1);
y = linspace(0,1,N(2)+1)';
y = y(1:end-1);
z = linspace(0,1,N(3)+1)';
z = z(1:end-1);
[X,Y,Z] = meshgrid(x,y,z);

pl = rho(:,:,:,end);
[~,pos] = max(pl(:));
[m1, m2, m3] = ind2sub(size(pl),pos);
pad = 8;
xx = x(m1-pad:m1+pad);
yy = y(m2-pad:m2+pad);
zz = z(m3-pad:m3+pad);
[XX,YY,ZZ] = meshgrid(xx,yy,zz);

writerObj = VideoWriter('rho3d.avi');
open(writerObj);
set(gcf,'Renderer','zbuffer');
% Create a sequence of frames and write each frame to the file
for i = 1:length(a)
pl = rho(:,:,:,i);
subplot(2,3,1)
h = slice(X,Y,Z, pl, x, y, z, 'cubic');
set(h, 'EdgeColor','none', 'FaceColor','interp')
alpha(.1)
subplot(2,3,2)
pl1 = pl;
plm = mean(pl1(:));
pl1(abs(pl1 - plm)/abs(plm) < 0.95) = nan;
hh = slice(X,Y,Z,pl1, x, y, z, 'cubic');
set(hh, 'EdgeColor','none', 'FaceColor','interp')
alpha(.3)
title(['rho @ t = ' num2str(t(i)) ', a = ' num2str(a(i))]);
subplot(2,3,3)
pl2 = pl(m1-pad:m1+pad,m2-pad:m2+pad,m3-pad:m3+pad);
% g = slice(pl, m1, m2, m3, 'cubic');
g = slice(XX,YY,ZZ,pl2, x(m1), y(m2), z(m3), 'cubic');
set(g, 'EdgeColor','none', 'FaceColor','interp')
alpha(0.5)
subplot(2,3,4)
surf(y,z,log10(squeeze(pl(m1,:,:)))); xlabel('y'); ylabel('z'); zlabel('log10(rho)');
shading interp
subplot(2,3,5)
surf(x,z,log10(squeeze(pl(:,m2,:)))); xlabel('x'); ylabel('z'); zlabel('log10(rho)');
shading interp
subplot(2,3,6)
surf(x,y,log10(squeeze(pl(:,:,m3)))); xlabel('x'); ylabel('y'); zlabel('log10(rho)');
shading interp

frame = getframe(gcf);
writeVideo(writerObj,frame);
end
close(writerObj);

system('ffmpeg -i rho3d.avi -vcodec libx264 -filter:v "setpts=15.0*PTS" rho3d2.mp4')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions from karsten

name = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/karsten/data_64.dat';
data = importdata(name);
data = squeeze(data(:,4:5));
phi = reshape(data(:,1),64,64,64);
dphi = reshape(data(:,2),64,64,64);
rho = dphi.^2 / 2 + mass^2 * phi.^2 / 2;

writerObj = VideoWriter('ic3d.avi');
open(writerObj);
set(gcf,'Renderer','zbuffer');
for i = 1:64
subplot(3,3,1)
surf(squeeze(phi(i,:,:)))
title('phi projected to yz')
subplot(3,3,2)
surf(squeeze(phi(:,i,:)))
title('phi projected to xz')
subplot(3,3,3)
surf(squeeze(phi(:,:,i)))
title('phi projected to xy')

subplot(3,3,4)
surf(squeeze(dphi(i,:,:)))
title('dphi projected to yz')
subplot(3,3,5)
surf(squeeze(dphi(:,i,:)))
title('dphi projected to xz')
subplot(3,3,6)
surf(squeeze(dphi(:,:,i)))
title('dphi projected to xy')

subplot(3,3,7)
surf(squeeze(rho(i,:,:)))
title('rho projected to yz')
subplot(3,3,8)
surf(squeeze(rho(:,i,:)))
title('rho projected to xz')
subplot(3,3,9)
surf(squeeze(rho(:,:,i)))
title('rho projected to xy')
frame = getframe(gcf);
writeVideo(writerObj,frame);
end
close(writerObj);