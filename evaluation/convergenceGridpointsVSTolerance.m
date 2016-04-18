method = 'hyp';
nums = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
gridpoints = [32, 64, 96, 128, 192, 256];
prefix= 'check3/';
suffix = '';

% name = [prefix method '_' num2str(max(gridpoints)) '_ref'];
% loadData
% rkphii = phi(:,1);
% rkphi = phi(:,end);
% rkx = x;

errsall = zeros(length(nums), length(gridpoints));
stepsall = zeros(size(errsall));

for ll = 1:length(gridpoints)

    Ngrid = gridpoints(ll);
    basename = [prefix method '_' num2str(Ngrid) '_'];

    name = [basename 'ref'];
    loadData
    rkphii = phi(:,1);
    rkphi = phi(:,end);
    rkx = x;
    
    errs1 = zeros(length(nums), 1);
    % errs2 = zeros(length(nums), 1);
    runtimes1 = zeros(length(nums), 1);
    % runtimes2 = zeros(length(nums), 1);
    steps1 = zeros(length(nums), 1);
    % steps2 = zeros(length(nums), 1);
    % phisi = {};
    % phisf = {};
    for kk = 1:length(nums)
        name = [basename num2str(nums(kk)) suffix];
        loadData
%         phitmp = spline(rkx, rkphi, x);
%         errs1(kk) = norm(phi(:,end) - phitmp) / N;
        errs1(kk) = norm(phi(:,end) - rkphi) / N;
%         runtimes1(kk) = runtime;
        steps1(kk) = steps;
    end
    errsall(:,ll) = errs1(:);
    stepsall(:,ll) = steps1(:);
end

figure
semilogy(nums, errsall, nums, 1e-3 * 10.^(-nums),'linewidth',2); xlabel('rel. tol. 10^{-x}'); ylabel('err');
title([method ': error wrt reference']);
legend('32','64','96','128','192','256');

figure
plot(nums, stepsall,'linewidth',2); xlabel('rel. tol. 10^{-x}'); ylabel('total steps');
title([method ': total number of steps']);
legend('32','64','96','128','192','256');

figure
semilogy(gridpoints, errsall,'linewidth',2); xlabel('gridpoints'); ylabel('err');
title([method ': error wrt reference']);
legend('3','4','5','6','7','8','9','10','11','12','13');