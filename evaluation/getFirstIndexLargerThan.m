function [idx, warn] = getFirstIndexLargerThan(f, val)
    warn = 0;
    idx = find(f(:) > val, 1);
    if isempty(idx)
        idx = length(f(:));
    end
    if idx == 1 || idx == length(f(:))
        disp 'Index at boundaries!'
        warn = 1;
    end
end