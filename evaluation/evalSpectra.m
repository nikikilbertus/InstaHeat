%% evaluate spectra of multiple runs

%% setup (user input)
% setup parameterization for the runs
ind = [512 1024 1536 2048 3072];
% ind = [0.001, 0.0005, 0.0001, 0.00005, 0.00001];
% ind = 2.^(7:16);
% construct the file names: prefix, suffix, indexset
pre = 'massres1d/';
suf = '_64_0.001';
% plot spectra at
as = [];
rhormss = [0.5, 1, 2, 3, 4, 5];
%figure offset
os = 0;

%% run

loadDsets;
require('rhoS','psiS','phiS','phips','dphiS','rhops','psips','N','spatial_bounds_x','mass');
nn = length(ind);
colors = jet;
ncol = size(colors,1);
colors = colors(1:int64(floor(ncol/nn)):end,:);
set(0,'DefaultAxesColorOrder',colors)
idx = zeros(length(as) + length(rhormss),2);
idxa =  zeros(nn, size(idx,1));
idxval = zeros(nn, 2, size(idx,1));
legendinfo = cell(nn,1);

for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;

    nbins = size(phips,2);
    N = N(1);
    L = spatial_bounds_x(2)-spatial_bounds_x(1);
    kmin = 2*pi/L; kmax = sqrt(3)*kmin*N/2;
    k = linspace(kmin,kmax,nbins);

    legendinfo{i} = num2str(ind(i));

    figure(os+1);
    subplot(1,2,1); loglog(a,rhorms); hold on;
    subplot(1,2,2); loglog(a,rhomean); hold on;

    figure(os+2);
    subplot(1,2,1); loglog(a,max(abs(psimin),abs(psimax))); hold on;
    subplot(1,2,2); loglog(a,psistd); hold on;

    figure(os+3);
    subplot(2,2,1); loglog(k,phips(2,:)); hold on;
    subplot(2,2,2); loglog(k,rhops(2,:)); hold on;
    subplot(2,2,3); loglog(k,phips(2,:)/phips(2,1)); hold on;
    subplot(2,2,4); loglog(k,rhops(2,:)/rhops(2,1)); hold on;
    
    for kk = 1:length(as)
        [r1,r2] = getIndexClosestTo(a,as(kk));
        idx(kk,:) = [r1 r2];
    end
    for kk = 1:length(rhormss)
        [r1,r2] = getFirstIndexLargerThan(rhorms,rhormss(kk));
        idx(length(as)+kk,:) = [r1 r2];
    end
    [~, perm] = sort(idx(:,1));
    idx = idx(perm,:);

    for kk = 1:size(idx,1)
        figure(os+3+kk);
        if idx(kk,2) ~= 1
            idxa(i,kk) = a(idx(kk,1));
            idxval(i,1,kk) = rhorms(idx(kk,1));
            idxval(i,2,kk) = rhomean(idx(kk,1));
            subplot(2,2,1); loglog(k,phips(idx(kk,1),:)); hold on;
            subplot(2,2,2); loglog(k,rhops(idx(kk,1),:)); hold on;
            subplot(2,2,3); loglog(k,phips(idx(kk,1),:)/phips(idx(kk,1),1)); hold on;
            subplot(2,2,4); loglog(k,rhops(idx(kk,1),:)/rhops(idx(kk,1),1)); hold on;
        else
            idxa(i,kk) = nan;
            idxval(i,:,kk) = nan;
            subplot(2,2,1); loglog(k,k*nan); hold on;
            subplot(2,2,2); loglog(k,k*nan); hold on;
            subplot(2,2,3); loglog(k,k*nan); hold on;
            subplot(2,2,4); loglog(k,k*nan); hold on;
        end
    end

    pot = mass^2*(phistd.^2+phimean.^2)/2;
    temp = (dphistd.^2+dphimean.^2)/2;
    grad = rhomean - pot - temp;
    resid = abs(rhomean - temp - grad - pot);
    figure(os+size(idx,1)+4)
    subplot(2,2,1); loglog(a,pot); hold on;
    subplot(2,2,2); loglog(a,temp); hold on;
    subplot(2,2,3); loglog(a,grad); hold on;
    subplot(2,2,4); loglog(a,resid); hold on;
    
    display(sprintf('processed %i of %i: %s', i, nn, name));
end

figure(os+1);
subplot(1,2,1);
for kk = 1:length(as)
    loglog([as(perm(kk)) as(perm(kk))],get(gca, 'ylim'),'b:');
end
for kk = 1:length(rhormss)
    loglog(get(gca, 'xlim'),[rhormss(perm(kk)) rhormss(perm(kk))],'b:');
end
for i = 1:nn
    loglog(idxa(i,:), squeeze(idxval(i,1,:)),'rx');
end
xlabel('a'); ylabel('std(\rho) / |<\rho>|'); legend(legendinfo,'location','northwest'); shg;
subplot(1,2,2);
for i = 1:nn
    loglog(idxa(i,:), squeeze(idxval(i,2,:)),'rx');
end
xlabel('a'); ylabel('<\rho>'); legend(legendinfo,'location','northeast'); shg;
figure(os+2);
subplot(1,2,1); xlabel('a'); ylabel('absmax \psi'); legend(legendinfo); shg;
subplot(1,2,2); xlabel('a'); ylabel('std(\psi)'); legend(legendinfo); shg;
figure(os+3);
subplot(2,2,1); xlabel('k'); ylabel('init powspec phi'); legend(legendinfo,'location','southwest'); shg;
subplot(2,2,2); xlabel('k'); ylabel('init powspec rho'); legend(legendinfo,'location','southwest'); shg;
for kk = 1:size(idx,1)
    figure(os+3+kk);
    subplot(2,2,1); xlabel('k'); ylabel('powspec phi'); legend(legendinfo,'location','southwest'); shg;
    title(['number ' num2str(kk)])
    subplot(2,2,2); xlabel('k'); ylabel('powspec rho'); legend(legendinfo,'location','southwest'); shg;
end
figure(os+size(idx,1)+4);
subplot(2,2,1); xlabel('a'); ylabel('<V>'); legend(legendinfo); shg;
subplot(2,2,2); xlabel('a'); ylabel('temporal E');
subplot(2,2,3); xlabel('a'); ylabel('gradient E');
subplot(2,2,4); xlabel('a'); ylabel('residual');