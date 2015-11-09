Nfiles = 4;
Nx = 501;
sol = @(x,t) zeros(size(x)); %sin(x) .* cos(t);
dtf = @(filenum) 0.1 * 2.^(-filenum);
finaltime = 10;

choice = questdlg('Show plots?', 'Display options', 'Yes', 'No', 'Yes');

switch choice
    case 'Yes'
         plotopt = 1;
    case 'No'
        plotopt = 0;
end

errs = zeros(Nfiles, 1);
energyloss = zeros(Nfiles, 1);
for i = 0:Nfiles-1
    [errs(i+1), energyloss(i+1)] = ...
        analyzeResult(i, Nx, sol, finaltime, dtf, plotopt);
end
dt = dtf(linspace(0, Nfiles-1, Nfiles));
dt = flip(dt);

% slope = logfit(dt', errs, 'loglog');
% title(['error on Tf depending on dt, exp = ' num2str(slope)])
% shg; pause();

slope = logfit(dt', energyloss, 'loglog');
title(['power = ' num2str(slope)])
xlabel('dt')
ylabel('energyloss')
shg; pause();