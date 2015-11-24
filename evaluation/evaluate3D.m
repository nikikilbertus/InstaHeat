Nfiles = 1;
dtf = @(filenum) 0.1 * 2.^(-filenum); %0.1 -filenum *0.01;
Nx = 32;
Ny = 32;
Nz = 32;
%dtf = @(filenum) 2 * pi / Nx / 10;
Ntot = Nx * Ny * Nz;
fieldWriteOutSize = 100;
powSpecWriteOutSize = 100;

fields = zeros(Nfiles, fieldWriteOutSize * Ntot);
for num = 0:Nfiles-1
    dt = dtf(num);
    prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';
    name = [prefix 'a_00' int2str(num) '.txt'];
    frwa = csvread(name);

    Nt = length(frwa);
    t = linspace(0, (Nt - 1) * dt, Nt);
    
    name = [prefix 'rho_00' int2str(num) '.txt'];
    rho = csvread(name);
    H = sqrt(rho / 3);
%     I = (t>0.1) & (t<0.5);
%     tfit = t(I);
%     loglog(t, rho, tfit, tfit.^(-1) * max(rho(I)) / tfit(1)^(-1));
    
    name = [prefix 'field_00' int2str(num) '.txt'];
    phi = csvread(name);
    phi = reshape(phi, Ntot, length(phi) / Ntot);
    phiavg = mean(phi);
    plot(phiavg);
    xlabel('time');
    ylabel('<phi>');
    shg;
    pause;
    
    if length(frwa) ~= length(rho)
       error('somethings wrong, check the parameters!') 
    end

%     [max,I] = findpeaks(rho);
%     loglog(frwa, rho);
%     hold on;
%     slopeup = logfit(frwa(I), rho(I), 'loglog');
%     slopeup = logfit(frwa, rho, 'loglog');
%     title(['dt = ' num2str(dt) ', slope = ' num2str(slopeup)]);
%     xlabel('a')
%     ylabel('rho');
%     shg;
%     pause;
    
    loglog(frwa, rho, frwa, frwa.^(-4) * rho(1));
    title(['dt = ' num2str(dt) ', slope = -4']);
    xlabel('a')
    ylabel('rho');
    shg;
    pause;
    
%     slopeup = logfit(frwa, H, 'loglog');
%     title(['dt = ' num2str(dt) ', slope = ' num2str(slopeup)]);
%     xlabel('a')
%     ylabel('H');
%     shg;
%     pause;
    
    loglog(frwa, H, frwa, frwa.^(-2) * H(1));
    title(['dt = ' num2str(dt) ', slope = -2']);
    xlabel('a')
    ylabel('H');
    shg;
    pause;
end

Nt = length(rho);
prefix = '~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/';
name = [prefix 'pow_spec_000.txt'];
powspec = csvread(name);
powspec = reshape(powspec, length(powspec) / (powSpecWriteOutSize), powSpecWriteOutSize);
for j=1:2
    surf(log(powspec(j:end, : )+1e-10));
%     set(gca, 'ZScale', 'log');
    shading interp; lighting phong;
    zlabel('log')
    ylabel('|k|');
    xlabel('nt');
    shg;
    pause;
    contourf
end

if powSpecWriteOutSize == fieldWriteOutSize
    parseval = zeros(1, fieldWriteOutSize);
    for i = 1:fieldWriteOutSize
    parseval(i) = abs(sqrt(sum(powspec(:,i))) - norm(phi(:, i)));
    end
    plot(parseval);
    title(['parseval, max error = ' num2str(max(parseval))]);
    xlabel('t');
    ylabel('||phi(k)|| - ||phi(x)||');
    shg;
    pause;
end

if Nfiles > 3
    dt = dtf(linspace(0, Nfiles-1, Nfiles));

    dnormsf = zeros(Nfiles - 1, 1);
    for i = 1:Nfiles-1
       dnormsf(i) = norm( fields(i+1, : ) - fields(i, : ) );
    end
    
    I = dnormsf > 0;
    I(1) = 0;
    dnormsfeff = dnormsf(I);
    ddt = zeros(5, Nfiles - 1);
    for i = 1:4
        ddt(i, : ) = abs(diff(dt.^i));
        ddteff = ddt(i, I);
        plot(ddteff, dnormsfeff, '-o');
        xlabel(['dt(i+1)^' num2str(i) ' - dt(i)^' num2str(i)]);
        ylabel('err(i+1) - err(i)');
        shg;
        pause();
    end
end