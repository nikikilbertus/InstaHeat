%% evaluate multiple runs with different resolutions or cutoffs

%% setup (user input)
% setup different gridpoints or cutoffs in an array
% res = [1 2 3 4];
% res = [0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001];
res = 2.^(7:16);
% construct the file names: prefix, suffix, indexset
pre = 'masses/1d/5e5/1024_64_';
suf = '';
ind = res;
%figure offset
os = 0;

%% run

loadDsets;
require(dsetsSummary,'constraints','steps_total','t','phips','rhops','psips','N','spatial_bounds_x');
nn = length(res);
colors = jet;
ncol = size(colors,1);
set(0,'DefaultAxesColorOrder',colors(1:int64(floor(ncol/nn)):end,:))
rhos = zeros(nn,2); % rhorms values at beginning and end
steps = zeros(nn,1);

for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;
    
    try
    steps(i) = steps_total;
    catch
    end
    
    nbins = size(phips,2);
    N = N(1);
    L = spatial_bounds_x(2)-spatial_bounds_x(1);
    kmin = 2*pi/L; kmax = sqrt(3)*kmin*N/2;
    k = (1:nbins)*kmax/nbins;
    
    rhos(i,:) = [rhorms(1); rhorms(end)];
    legendinfo{i} = ['N=' num2str(res(i))];
    
    figure(os+1); plot(diff(t)); hold on;
    
    figure(os+2); loglog(a,rhorms); hold on;
    
    figure(os+3);
    subplot(1,2,1); loglog(a,abs(phimean)); hold on;
    subplot(1,2,2); loglog(a,phistd); hold on;
    
    figure(os+4);
    subplot(1,2,1); loglog(a,abs(dphimean)); hold on;
    subplot(1,2,2); loglog(a,dphistd); hold on;
    
    figure(os+5);
    subplot(1,2,1); loglog(a,max(abs(psimin),abs(psimax))); hold on;
    subplot(1,2,2); loglog(a,psistd); hold on;
    
    figure(os+6);
    subplot(1,2,1); loglog(a,max(abs(dpsimin),abs(dpsimax))); hold on;
    subplot(1,2,2); loglog(a,dpsistd); hold on;
    
    figure(os+7);
    subplot(1,2,1); loglog(a,rhomean); hold on;
    subplot(1,2,2); loglog(t,a); hold on;
    
    if exist('constraints','var')
        figure(os+8);
        subplot(1,2,1); loglog(a, abs(constraints(:,1))); hold on;
        subplot(1,2,2); loglog(a, abs(constraints(:,2))); hold on;
        figure(os+9);
        subplot(1,2,1); loglog(a, abs(constraints(:,3))); hold on;
        subplot(1,2,2); loglog(a, abs(constraints(:,4))); hold on;
    end
    
    figure(os+10);
    subplot(1,2,1); loglog(k,phips(2,:)); hold on;
    subplot(1,2,2); loglog(k,rhops(1,:)); hold on;
    
    figure(os+11);
    [idx, warn] = getIndexClosestTo(rhorms, 1);
    if warn ~= 1
        subplot(1,2,1); loglog(k,phips(idx,:)); hold on;
        subplot(1,2,2); loglog(k,rhops(idx,:)); hold on;
    else
        subplot(1,2,1); loglog(k,1); hold on;
        subplot(1,2,2); loglog(k,1); hold on;
    end
    
    display(sprintf('processed %i of %i: %s', i, nn, name));
end

figure(os+1); xlabel('#step'); ylabel('\Delta t'); legend(legendinfo,'location','southeast'); shg;
figure(os+2);
loglog([a(1) a(end)], [1 1],':');
xlabel('a'); ylabel('std(\rho) / |<\rho>|'); legend(legendinfo,'location','northwest'); shg;
figure(os+3);
subplot(1,2,1); xlabel('a'); ylabel('|<\phi>|'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(\phi)'); legend(legendinfo); shg;
figure(os+4);
subplot(1,2,1); xlabel('a'); ylabel('|<d\phi>|'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(d\phi)'); legend(legendinfo); shg;
figure(os+5);
subplot(1,2,1); xlabel('a'); ylabel('absmax \psi'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(\psi)'); legend(legendinfo); shg;
figure(os+6);
subplot(1,2,1); xlabel('a'); ylabel('absmax d\psi'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std d\psi'); legend(legendinfo); shg;
figure(os+7);
subplot(1,2,1); xlabel('a'); ylabel('rhomean'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('t'); ylabel('a'); legend(legendinfo,'location','northwest'); shg;
if exist('constraints','var')
    figure(os+8);
    subplot(1,2,1); xlabel('a'); ylabel('ham cstr l2'); legend(legendinfo); shg;
    subplot(1,2,2); xlabel('a'); ylabel('ham cstr \infty'); legend(legendinfo); shg;
    figure(os+9);
    subplot(1,2,1); xlabel('a'); ylabel('mom cstr l2'); legend(legendinfo); shg;
    subplot(1,2,2); xlabel('a'); ylabel('mom cstr \infty'); legend(legendinfo); shg;
end
figure(os+10);
subplot(1,2,1); xlabel('k'); ylabel('init powspec phi'); legend(legendinfo,'location','southwest'); shg;
subplot(1,2,2); xlabel('k'); ylabel('init powspec rho'); legend(legendinfo,'location','southwest'); shg;
figure(os+11);
subplot(1,2,1); xlabel('k'); ylabel('final powspec phi'); legend(legendinfo,'location','southwest'); shg;
subplot(1,2,2); xlabel('k'); ylabel('final powspec rho'); legend(legendinfo,'location','southwest'); shg;
figure
loglog(res.^3, rhos, '-o'); xlabel('N^3'); ylabel('std(\rho) / |<\rho>|'); shg;
legend('initial','final');
figure
loglog(res.^3, steps, '-o'); xlabel('N^3'); ylabel('steps'); shg;
set(0,'DefaultAxesColorOrder','remove')
