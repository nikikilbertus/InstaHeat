%% evaluate runs with different masses

%% setup parameters
% set up 'masses' array
ml = 3; mh = 4;
masses = 5*10.^(-(ml:mh));
% construct the file names: prefix, suffix, indexset
pre = 'gw4/64_5e-';
suf = '_2e4';
ind = ml:mh;

%% run
require('rhorms');
nn = mh-ml+1;
rhormsi = zeros(nn,4);
maxrhos = zeros(nn,1);
for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;
    maxrhos(i) = max(rhorms);
    loglog(a,rhorms/rhorms(1),'linewidth',2); hold on
    legendinfo{i} = ['mass=' num2str(masses(i))];
    rhormsi(i,1) = rhorms(1);
    [~,idx] = min((a - a(end)/3).^2);
    rhormsi(i,2) = rhorms(idx);
    [~,idx] = min((a - 2*a(end)/3).^2);
    rhormsi(i,3) = rhorms(idx);
    rhormsi(i,4) = rhorms(end);
end
hold off; xlabel('a'); ylabel('std(\rho)/|<\rho>|');
legend(legendinfo,'location','southeast');
figure
loglog(masses,rhormsi, masses,masses/masses(1)*rhormsi(1,1)*0.9,'--','linewidth',2);
xlabel('mass'); ylabel('std(\rho)/|<\rho>|');
legend('a=1', ['a=' num2str(a(end)/3)], ['a=' num2str(2*a(end)/3)], ...
    ['a=' num2str(a(end))],'linear reference','location','northwest');