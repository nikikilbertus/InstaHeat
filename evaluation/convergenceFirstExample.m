basename = 'resolution_rk_6_amp';
nums = [4, 5, 6, 7, 8, 9];

as = {};
dpsis = {};
dpsifd2s = {};
for kk = 1:length(nums)
    name = [basename num2str(nums(kk))];
    loadData
    as{kk} = a;
    dpsis{kk} = dpsi;
    dpsifd2s{kk} = dpsifd2;
end

Nh = round(N/2);
figure
for kk = 1:length(nums)
    tmp = dpsis{kk};
    plot(as{kk}, log10(tmp(Nh,:)))
    hold on
end
hold off
xlabel('a')
ylabel('dpsi middle')
legend('1e-4', '1e-5','1e-6','1e-7','1e-8','1e-9')

figure
for kk = 1:length(nums)
    dpsitmp = dpsis{kk};
    dpsitmpfd = dpsifd2s{kk};
    atmp = as{kk};
    plot(atmp(2:end-1), log10(abs(dpsitmp(Nh,2:end-1) - dpsitmpfd(Nh,:))))
    hold on
end
hold off
xlabel('a')
ylabel('difference dpsi')
legend('1e-4', '1e-5','1e-6','1e-7','1e-8','1e-9')

figure
for kk = 1:length(nums)
    dpsitmp = dpsis{kk};
    dpsitmpfd = dpsifd2s{kk};
    atmp = as{kk};
    plot(atmp(2:end-60000), log10(abs(dpsitmp(Nh,2:end-60000) - dpsitmpfd(Nh,1:end-59999))))
    hold on
end
hold off
xlabel('a')
ylabel('difference dpsi')
legend('1e-4', '1e-5','1e-6','1e-7','1e-8','1e-9')

for kk = 1:length(nums)
    dpsitmp = dpsis{kk};
    dpsitmpfd = dpsifd2s{kk};
    atmp = as{kk};
    figure
    surf(atmp(2:end-60000),x,dpsitmp(:,2:end-60000) - dpsitmpfd(:,1:end-59999))
    shading interp
    title(['diff dpsi 1e-' num2str(kk + 3)])
end