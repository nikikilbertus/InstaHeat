Nfiles = 4;
Nx = 101;
sol = @(x,t) zeros(size(x)); %sin(x) .* cos(t);
dtf = @(filenum) 0.01 * 2.^(-filenum); %0.1 -filenum *0.01;
finaltime = 1.5;

[~, ~, x] = mkSpectralOperators({Nx, 'fourier'});

choice = questdlg('Show plots?', 'Display options', 'Yes', 'No', 'Yes');

switch choice
    case 'Yes'
         plotopt = 1;
    case 'No'
        plotopt = 0;
end

errs = zeros(Nfiles, 1);
energyloss = zeros(Nfiles, 1);
Tf = zeros(Nfiles, Nx);
for i = 0:Nfiles-1
    [errs(i+1), energyloss(i+1), Tf(i+1, : )] = ...
        analyzeResult(i, Nx, x, sol, finaltime, dtf, plotopt);
end
dt = dtf(linspace(0, Nfiles-1, Nfiles));

% slope = logfit(dt', errs, 'loglog');
% title(['error on Tf depending on dt, exp = ' num2str(slope)])
% shg; pause();

slope = logfit(dt', energyloss, 'loglog');
title(['power = ' num2str(slope)])
xlabel('dt')
ylabel('energyloss')
shg; pause();

dnormsf = zeros(Nfiles-1,1);
for i = 1:Nfiles-1
   dnormsf(i) = norm( Tf(i+1, : ) - Tf(i, : ) );
end

ddt = zeros(5,Nfiles-1);
for i = 1:5
    ddt(i, : ) = abs(diff(dt.^i));
    plot(ddt(i,:), dnormsf, 'o');
    xlabel(['dt(i+1)^' num2str(i) ' - dt(i)^' num2str(i)]);
    ylabel('err(i+1) - err(i)');
    shg;
    pause();
    plot(dnormsf'./ddt(i,:), 'o');
    ylabel(['(err(i+1) - err(i)) / (dt(i+1)^' num2str(i) ' - dt(i)^' num2str(i) ')']);
    shg;
    pause();
end