method = 'par';
nums = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
gridpoints = [32, 64, 96, 128];
suffix = '';

errsall = zeros(length(nums), length(gridpoints));
stepsall = zeros(size(errsall));
for ll = 1:length(gridpoints)

    Ngrid = gridpoints(ll);
    name = ['check/' method '_' num2str(Ngrid) '_ref'];
    loadData
    rkphii = phi(:,1);
    rkphi = phi(:,end);

    basename = ['check/' method '_' num2str(Ngrid) '_'];

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
    %     phisi{kk} = phi(:,1);
    %     phisf{kk} = phi(:,end);
        errs1(kk) = norm(phi(:,end) - rkphi(:)) / N;
%         runtimes1(kk) = runtime;
        steps1(kk) = steps;
    %     name = [basename num2str(nums(kk)) '_a0'];
    %     loadData
    %     errs2(kk) = norm(phi(:,end) - rkphi(:)) / N;
    %     runtimes2(kk) = runtime;
    %     steps2(kk) = steps;
    end
    errsall(:,ll) = errs1(:);
    stepsall(:,ll) = steps1(:);
end

figure
semilogy(nums, errsall,'linewidth',2); xlabel('rel. tol. 10^{-x}'); ylabel('err');
title('error wrt reference');
legend('32','64','96','128');

figure
plot(nums, stepsall,'linewidth',2); xlabel('rel. tol. 10^{-x}'); ylabel('total steps');
title('total number of steps');
legend('32','64','96','128');

figure
semilogy(gridpoints, errsall,'linewidth',2); xlabel('gridpoints'); ylabel('err');
title('error wrt reference');
legend('3','4','5','6','7','8','9','10','11','12','13');