function [ps] = mkPowerSpectrum(f, bins)

Nall = numel(f);
fmom = fftn(f);
fmom = abs(fmom).^2;

x = arrayfun(@(in) colon(1,in), size(f), 'UniformOutput',false);
if ndims(f) == 1
    kks = x{1}.^2;
elseif ndims(f) == 2
    [X, Y] = meshgrid(x{1}, x{2});
    kks = X.^2 + Y.^2;
elseif ndims(f) == 3
    [X, Y, Z] = meshgrid(x{1}, x{2}, x{3});
    kks = X.^2 + Y.^2 + Z.^2;
end

[~, I] = sort(kks(:));
fmom = fmom(I);

perbin = ceil(Nall / bins)
newperbin = perbin;
ps = zeros(1, bins);
inds = zeros(1,Nall);
for i = 1:bins
    tmp = 0;
    if i == bins && bins*perbin > Nall
        newperbin = bins * perbin - Nall;
        disp('perbin set to')
        disp(newperbin)
    end
    for j = 1:min(perbin, newperbin)
        tmp = tmp + fmom((i-1)*perbin + j);
        inds((i-1)*perbin + j) = (i-1)*perbin + j;
    end
    ps(i) = tmp;
end


plot(inds)
shg
figure
return

kkmax = max(kks(:));
dkk = kkmax / bins;
% ps = zeros(1, kkmax);

for i = 1:Nall
%     ps(kks(i)) = fmom(i)/Nall;
    idx = int64(ceil(kks(i) / dkk));
    if(idx > bins)
        error('wrong index')
    end
    ps(idx) = ps(idx) + fmom(i)/Nall;
end

end