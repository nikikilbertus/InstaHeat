%% evaluate runs with different masses

%% setup parameters
% set up 'masses' array
% ml = 3; mh = 4;
masses = [0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001];
% construct the file names: prefix, suffix, indexset
pre = 'masses/48_32_';
suf = '';
ind = masses;

%% run
require('rhoS');
nn = length(masses);
rhormsi = zeros(nn,4);
maxrhos = zeros(nn,1);
for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;
    maxrhos(i) = max(rhorms);
    subplot(1,2,1)
    loglog(a,rhorms,'linewidth',2); hold on
    subplot(1,2,2)
    loglog(a,rhorms/rhorms(1),'linewidth',2); hold on
    legendinfo{i} = ['mass=' num2str(masses(i))];
    rhormsi(i,1) = rhorms(1); % rhorms in the beginning
    [~,idx] = min((a - a(end)/10).^2);
    rhormsi(i,2) = rhorms(idx); % rhorms at one tenth of the evolution
    [~,idx] = min((a - a(end)/2).^2);
    rhormsi(i,3) = rhorms(idx); % rhorms at half of the evolution
    rhormsi(i,4) = rhorms(end); % rhorms in the end
    display(sprintf('processed %i of %i', i, nn));
end
hold off; 
subplot(1,2,1)
xlabel('a'); ylabel('std(\rho)/|<\rho>|');
legend(legendinfo,'location','southeast');
subplot(1,2,2)
xlabel('a'); ylabel('std(\rho)/|<\rho>| normalized');
legend(legendinfo,'location','southeast');
figure
loglog(masses,rhormsi, masses,masses/masses(1)*rhormsi(1,1)*0.9,'--','linewidth',2);
xlabel('mass'); ylabel('std(\rho)/|<\rho>|');
legend('a=1', ['a=' num2str(a(end)/3)], ['a=' num2str(2*a(end)/3)], ...
    ['a=' num2str(a(end))],'linear reference','location','northwest');