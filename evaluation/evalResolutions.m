%% evaluate multiple runs with different resolutions or cutoffs

%% setup (user input)
% setup different gridpoints or cutoffs in an array
res = [16 32 48 64 96];
% construct the file names: prefix, suffix, indexset
pre = 'resolution1e-3/';
suf = '_16_1e-3';
ind = res;

%% run
loadDsets;
require(dsetsSummary,'steps_total','t');
nn = length(res);
rhos = zeros(nn,2); % rhorms values at beginning and end
steps = zeros(nn,1);
for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;
    steps(i) = steps_total;
    rhos(i,:) = [rhorms(1); rhorms(end)];
    legendinfo{i} = ['N=' num2str(res(i))];
    figure(1); plot(diff(t)); hold on;
    figure(2); loglog(a,rhorms); hold on;
    figure(3);
    subplot(1,2,1); loglog(a,abs(phimean)); hold on;
    subplot(1,2,2); loglog(a,phistd); hold on;
    figure(4);
    subplot(1,2,1); loglog(a,abs(dphimean)); hold on;
    subplot(1,2,2); loglog(a,dphistd); hold on;
    figure(5);
    subplot(1,2,1); loglog(a,max(abs(psimin),abs(psimax))); hold on;
    subplot(1,2,2); loglog(a,psistd); hold on;
    figure(6);
    subplot(1,2,1); loglog(a,max(abs(dpsimin),abs(dpsimax))); hold on;
    subplot(1,2,2); loglog(a,dpsistd); hold on;
    
    if exist('hamcstrl2','var')
        figure(7); loglog(a, abs(hamcstrl2)/res(i)^3); hold on;
    end
    if exist('hamcstrinf','var')
        figure(8); loglog(a, abs(hamcstrinf)); hold on;
    end
    display(sprintf('processed %i of %i', i, nn));
end
for i = 1:6
    figure(i); hold off;
end
figure(1); xlabel('#step'); ylabel('\Delta t'); legend(legendinfo); shg;
figure(2); xlabel('a'); ylabel('std(\rho) / |<\rho>|'); legend(legendinfo); shg;
figure(3);
subplot(1,2,1); xlabel('a'); ylabel('|<\phi>|'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(\phi)'); legend(legendinfo); shg;
figure(4);
subplot(1,2,1); xlabel('a'); ylabel('|<d\phi>|'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(d\phi)'); legend(legendinfo); shg;
figure(5);
subplot(1,2,1); xlabel('a'); ylabel('absmax \psi'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(\psi)'); legend(legendinfo); shg;
figure(6);
subplot(1,2,1); xlabel('a'); ylabel('absmax d\psi'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std d\psi'); legend(legendinfo); shg;
if exist('hamcstrl2','var')
    figure(7); hold off; xlabel('a'); ylabel('ham cstr l2'); legend(legendinfo); shg;
end
if exist('hamcstrinf','var')
    figure(8); hold off; xlabel('a'); ylabel('ham cstr \infty'); legend(legendinfo); shg;
end
figure
loglog(res.^3, rhos, '-o'); xlabel('N^3'); ylabel('std(\rho) / |<\rho>|'); shg;
legend('initial','final');
figure
loglog(res.^3, steps, '-o'); xlabel('N^3'); ylabel('steps'); shg;