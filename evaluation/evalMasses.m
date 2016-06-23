%% evaluate multiple runs with different masses

%% setup (user input)
% set up 'masses' array
masses = [0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001];
% if those are long runs: how many time slices do we want (<0 for all)
nout = 3000;
% construct the input file names: prefix, suffix, indexset
pre = 'masses/48_32_';
suf = '';
ind = masses;
% construct the output file names
preout = '48_32_';
sufout = '';
indout = masses;

%% run
nn = length(masses);
name = [pre num2str(ind(1)) suf];
readDsets;
as = [1, max(a)/10, max(a)/2, max(a)];
rhormsi = zeros(nn,4);
psistdi = zeros(nn,4);
for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    if nout > 0
        clear dsets
        clear readI
        readDsets;
        readI = getLogIndices(length(a), nout);
    end
    require('rhoS','psiS');
    readDsets;
    legendinfo{i} = ['mass=' num2str(masses(i))];
    figure(1);
    subplot(1,2,1); loglog(a,rhorms); hold on;
    subplot(1,2,2); loglog(a,rhorms/rhorms(1)); hold on;
    figure(2);
    subplot(1,2,1); loglog(a,abs(psistd)); hold on;
    subplot(1,2,2); loglog(a,abs(psistd/psistd(1))); hold on;
    rhormsi(i,1) = rhorms(1); % rhorms in the beginning
    psistdi(i,1) = psistd(1); % psistd in the beginning
    [~,idx] = min((a - as(2)).^2);
    rhormsi(i,2) = rhorms(idx); % rhorms at one tenth of the evolution
    psistdi(i,2) = psistd(idx); % psistd at one tenth of the evolution
    [~,idx] = min((a - a(3)).^2);
    rhormsi(i,3) = rhorms(idx); % rhorms at half of the evolution
    psistdi(i,3) = psistd(idx); % psistd at half of the evolution
    rhormsi(i,4) = rhorms(end); % rhorms in the end
    psistdi(i,4) = psistd(end); % psistd in the end
    display(sprintf('processed %i of %i', i, nn));
    % output
    nameout = [preout num2str(indout(i)) sufout '.csv'];
    psiabsmax = max(abs(psimin),abs(psimax));
    psiabsstd = abs(psistd);
    T = table(a,rhomean,rhorms,psiabsstd,psiabsmax);
    writetable(T,nameout);
end
figure(1);
subplot(1,2,1);
for i = 1:4
    loglog([as(i) as(i)], [1e-2 5e0],'black');
end
xlabel('a'); ylabel('std(\rho) / |<\rho>|'); legend(legendinfo,'location','southeast');
subplot(1,2,2);
xlabel('a'); ylabel('std(\rho) / |<\rho>| normalized'); legend(legendinfo,'location','southeast');
figure(2);
subplot(1,2,1);
for i = 1:4
    loglog([as(i) as(i)], [1e-3 1e-1],'black');
end
xlabel('a'); ylabel('std(\psi)'); legend(legendinfo,'location','southeast');
subplot(1,2,2);
xlabel('a'); ylabel('std(\psi) normalized'); legend(legendinfo,'location','southeast');
% various intermediate rohrms values as function of masses
figure
loglog(masses,rhormsi, masses,masses/masses(1)*rhormsi(1,1)*0.9,'--','linewidth',2);
xlabel('mass'); ylabel('std(\rho) / |<\rho>|');
legend('a=1', ['a=' num2str(a(end)/3)], ['a=' num2str(2*a(end)/3)], ...
    ['a=' num2str(a(end))],'linear reference','location','northwest');
% various intermediate phistd values as function of masses
figure
loglog(masses,psistdi);
xlabel('mass'); ylabel('std(\psi)');
legend('a=1', ['a=' num2str(a(end)/3)], ['a=' num2str(2*a(end)/3)], ...
    ['a=' num2str(a(end))],'location','northwest');
return

%% output