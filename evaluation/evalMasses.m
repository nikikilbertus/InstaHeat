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
as = [1, max(a)/100, max(a)/10, max(a)];
rhormsi = zeros(nn,4);
rhoscali = zeros(nn,4);
psistdi = zeros(nn,4);
psimaxi = zeros(nn,4);
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
    psiabsmax = max(abs(psimin),abs(psimax));
    psiabsstd = abs(psistd);
    rhoscal = max(abs(rhomin),abs(rhomax))./rhomean;
    figure(1);
    subplot(1,3,1); loglog(a,rhorms); hold on;
    subplot(1,3,2); loglog(a,rhorms/rhorms(1)); hold on;
    subplot(1,3,3); loglog(a,rhoscal); hold on;
    figure(2);
    subplot(1,2,1); loglog(a,psiabsstd,a,psiabsmax); hold on;
    subplot(1,2,2); loglog(a,psiabsstd/psiabsstd(1)); hold on;
     % in the beginning
    rhormsi(i,1) = rhorms(1);
    rhoscali(i,1) = rhoscal(1);
    psistdi(i,1) = psistd(1);
    psimaxi(i,1) = psiabsmax(1);
    % at first intermediate point of the evolution
    [~,idx] = min((a - as(2)).^2);
    rhormsi(i,2) = rhorms(idx);
    rhoscali(i,2) = rhoscal(idx);
    psistdi(i,2) = psistd(idx);
    psimaxi(i,2) = psiabsmax(idx);
    % at second intermediate point of the evolution
    [~,idx] = min((a - as(3)).^2);
    rhormsi(i,3) = rhorms(idx);
    rhoscali(i,3) = rhoscal(idx);
    psistdi(i,3) = psistd(idx);
    psimaxi(i,3) = psiabsmax(idx); 
    % at end
    rhormsi(i,4) = rhorms(end);
    rhoscali(i,4) = rhoscal(end);
    psistdi(i,4) = psistd(end);
    psimaxi(i,4) = psiabsmax(end);
    display(sprintf('processed %i of %i', i, nn));
    % output
    nameout = [preout num2str(indout(i)) sufout '.csv'];
    T = table(a,rhomean,rhorms,rhoscal,psiabsstd,psiabsmax);
    writetable(T,nameout);
end
figure(1);
subplot(1,3,1);
for i = 1:4
    loglog([as(i) as(i)], [1e-3 5e0],'--black');
end
xlabel('a'); ylabel('std(\rho) / |<\rho>|'); legend(legendinfo,'location','southeast');
subplot(1,3,2);
xlabel('a'); ylabel('std(\rho) / |<\rho>| normalized'); legend(legendinfo,'location','southeast');
subplot(1,3,3);
xlabel('a'); ylabel('max(|\rho|)'); legend(legendinfo,'location','southeast');
figure(2);
subplot(1,2,1);
for i = 1:4
    loglog([as(i) as(i)], [1e-4 1e-2],'--black');
end
xlabel('a'); ylabel('std(\psi)'); legend(legendinfo,'location','southeast');
subplot(1,2,2);
xlabel('a'); ylabel('std(\psi) normalized'); legend(legendinfo,'location','southeast');
% various intermediate rohrms values as function of masses
figure
loglog(masses,rhormsi, masses,masses/masses(1)*rhormsi(1,1)*0.9,'--','linewidth',2);
xlabel('mass'); ylabel('std(\rho) / |<\rho>|');
legend('a=1', ['a=' num2str(as(2))], ['a=' num2str(as(3))], ...
    ['a=' num2str(a(end))],'linear reference','location','northwest');
% various intermediate phistd values as function of masses
figure
loglog(masses,psimaxi);
xlabel('mass'); ylabel('max(|\psi|)');
legend('a=1', ['a=' num2str(as(2))], ['a=' num2str(as(3))], ...
    ['a=' num2str(a(end))],'location','northwest');
return

T = table(masses',rhormsi,rhoscali,psistdi,psimaxi);
writetable(T,'intermediate.csv');