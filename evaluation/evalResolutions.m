%% evaluate multiple runs with different resolutions or cutoffs

%% setup (user input)
% setup different gridpoints or cutoffs in an array
res = [32 40 48 56 64 80 96 112 128 144 160 192 224 256];
% res = [32 40 48 56 64 80 96];
% construct the file names: prefix, suffix, indexset
pre = 'res2daggr/';
suf = '_1e-4_8';
ind = res;
%figure offset
os = 0;

%% run

loadDsets;
require(dsetsSummary,'constraints','steps_total','t');
nn = length(res);
colors = jet;
ncol = size(colors,1);
set(0,'DefaultAxesColorOrder',colors(1:int64(floor(ncol/nn)):end,:))
rhos = zeros(nn,2); % rhorms values at beginning and end
steps = zeros(nn,1);
for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;
    steps(i) = steps_total;
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
    display(sprintf('processed %i of %i', i, nn));
end
for i = (1:7)+os
    figure(os+i); hold off;
end
figure(os+1); xlabel('#step'); ylabel('\Delta t'); legend(legendinfo,'location','southeast'); shg;
figure(os+2); xlabel('a'); ylabel('std(\rho) / |<\rho>|'); legend(legendinfo,'location','northwest'); shg;
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
    subplot(1,2,1); xlabel('a'); ylabel('ham cstr l2'); legend(legendinfo); shg; hold off;
    subplot(1,2,2); xlabel('a'); ylabel('ham cstr \infty'); legend(legendinfo); shg; hold off;
    figure(os+9);
    subplot(1,2,1); xlabel('a'); ylabel('mom cstr l2'); legend(legendinfo); shg; hold off;
    subplot(1,2,2); xlabel('a'); ylabel('mom cstr \infty'); legend(legendinfo); shg; hold off;
end
figure
loglog(res.^3, rhos, '-o'); xlabel('N^3'); ylabel('std(\rho) / |<\rho>|'); shg;
legend('initial','final');
figure
loglog(res.^3, steps, '-o'); xlabel('N^3'); ylabel('steps'); shg;
set(0,'DefaultAxesColorOrder','remove')
