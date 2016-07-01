%% test how the root mean square of the same function differs on different grids
N = 16:8:96;
nn = length(N);
f = zeros(nn,1);
for i = 1:nn
    x = linspace(0,2*pi,N(i));
    [X,Y,Z] = meshgrid(x,x,x);
    tmp = sin(X) .* cos(Y) .* sin(2*Z) + 10;
    f(i) = std(tmp(:))/mean(tmp(:));
end
plot(f); shg;
(max(f) - min(f)) / min(f)