gridpoints = 128;

name = ['check/ell_' num2str(gridpoints) '_ref'];
loadData
rkphii = phi(:,1);
rkphi = phi(:,end);

basename = ['check/ell_' num2str(gridpoints) '_'];
nums = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];

errs1 = zeros(length(nums), 1);
% errs2 = zeros(length(nums), 1);
runtimes1 = zeros(length(nums), 1);
% runtimes2 = zeros(length(nums), 1);
steps1 = zeros(length(nums), 1);
% steps2 = zeros(length(nums), 1);
% phisi = {};
% phisf = {};
for kk = 1:length(nums)
    name = [basename num2str(nums(kk))];
    loadData
%     phisi{kk} = phi(:,1);
%     phisf{kk} = phi(:,end);
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
plot(nums, runtimes1); xlabel('rel. tol. 10^{-x}'); ylabel('runtime');
title('runtime of stepper');
legend('a=1e-14');

figure
plot(nums, steps1); xlabel('rel. tol. 10^{-x}'); ylabel('total steps');
title('total number of steps');
legend('a=1e-14');