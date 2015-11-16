Nfiles = 1;
dtf = @(filenum) 0.01 * 2.^(-filenum); %0.1 -filenum *0.01;

for num = 0:Nfiles-1
    dt = dtf(num);
    prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';
    name = [prefix 'a_00' int2str(num) '.txt'];
    frwa = csvread(name);

    name = [prefix 'energy_00' int2str(num) '.txt'];
    rho = csvread(name);

    if length(frwa) ~= length(rho)
       error('somethings wrong, check the parameters!') 
    end

%     [max,I] = findpeaks(rho);
%     loglog(frwa, rho);
%     hold on;
%     slopeup = logfit(frwa(I), rho(I), 'loglog');
    slopeup = logfit(frwa, rho, 'loglog');
    title(['dt = ' num2str(dt) ', slope = ' num2str(slopeup)]);
    xlabel('a')
    ylabel('rho');
%     legend('data', ['slope = ' num2str(slopeup)]);
%     hold off;
    shg;
    pause();
end