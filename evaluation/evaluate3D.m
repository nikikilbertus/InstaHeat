Nfiles = 1;
dtf = @(filenum) 0.01 * 2.^(-filenum); %0.1 -filenum *0.01;
Nx = 128;
Ny = 128;
Nz = 128;
Ntot = Nx * Ny * Nz;

fields = zeros(Nfiles, Ntot);
for num = 0:Nfiles-1
    dt = dtf(num);
    prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';
    name = [prefix 'a_00' int2str(num) '.txt'];
    frwa = csvread(name);

    name = [prefix 'rho_00' int2str(num) '.txt'];
    rho = csvread(name);

    name = [prefix 'field_00' int2str(num) '.txt'];
    fields(num + 1, : ) = csvread(name);
    
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

if Nfiles > 3
    dt = dtf(linspace(0, Nfiles-1, Nfiles));

    dnormsf = zeros(Nfiles - 1, 1);
    for i = 1:Nfiles-1
       dnormsf(i) = norm( fields(i+1, : ) - fields(i, : ) );
    end

    ddt = zeros(5, Nfiles - 1);
    for i = 1:4
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
end