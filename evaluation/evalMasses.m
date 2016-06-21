%% evaluate runs with different masses

%% setup parameters
% TODO
ml = 3; m = 4;
% masses = 1./(2*10.^(-(ml:m)));
masses = 5*10.^(-(ml:m));
rhormsi = zeros(m-ml+1,4);
maxrhos = zeros(m-ml+1,1);
for i = 1:m-ml+1
    name = ['gw4/64_5e-' num2str(i+ml-1) '_2e4'];
    evaluate3D
    maxrhos(i) = max(rhorms);
    loglog(a,rhorms/rhorms(1),'linewidth',2); hold on
    legendinfo{i} = ['mass=' num2str(masses(i))];
    rhormsi(i,1) = rhorms(1);
    [~,idx] = min((a - 1e2).^2);
    rhormsi(i,2) = rhorms(idx);
    [~,idx] = min((a - 1e3).^2);
    rhormsi(i,3) = rhorms(idx);
    rhormsi(i,4) = rhorms(end);
end
hold off; xlabel('a'); ylabel('std(\rho) / |<\rho>|');
legend(legendinfo,'location','southeast');
figure
loglog(masses, rhormsi, masses, masses/masses(1) * rhormsi(1,1) * 0.9,'--','linewidth',2);
xlabel('mass'); ylabel('std(\rho) / |<\rho>|');
legend('a=1', 'a=1e2', 'a=1e3', ['a=' num2str(a(end))],'linear reference','location','northwest');