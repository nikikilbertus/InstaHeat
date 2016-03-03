gridpoints = 192;

name = ['rk_' num2str(gridpoints)];
loadData
rkphi = phi(:,end);

basename = ['dp_' num2str(gridpoints) '_'];
nums = [3, 4, 5, 6, 7, 8, 9];

errs1 = zeros(length(nums), 1);
errs2 = zeros(length(nums), 1);
runtimes1 = zeros(length(nums), 1);
runtimes2 = zeros(length(nums), 1);
steps1 = zeros(length(nums), 1);
steps2 = zeros(length(nums), 1);
% phis = {};
for kk = 1:length(nums)
    name = [basename num2str(nums(kk))];
    loadData
%     phis{kk} = phi(:,end);
    errs1(kk) = norm(phi(:,end) - rkphi(:)) / N;
    runtimes1(kk) = runtime;
    steps1(kk) = steps;
%     name = [basename num2str(nums(kk)) '_a0'];
%     loadData
%     errs2(kk) = norm(phi(:,end) - rkphi(:)) / N;
%     runtimes2(kk) = runtime;
%     steps2(kk) = steps;
end

semilogy(nums, errs1); xlabel('rel. tol. 10^{-x}'); ylabel('err');
title('error wrt reference');

figure
plot(nums, runtimes1, nums, runtimes2); xlabel('rel. tol. 10^{-x}'); ylabel('runtime');
title('runtime of stepper');
legend('a=1e-14', 'a=0');

figure
plot(nums, steps1, nums, steps2); xlabel('rel. tol. 10^{-x}'); ylabel('total steps');
title('total number of steps');
legend('a=1e-14', 'a=0');