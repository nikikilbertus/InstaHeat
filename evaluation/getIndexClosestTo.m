function [idx, warn] = getIndexClosestTo(f, val)
    warn = 0;
    [~, idx] = min(abs(f(:) - val));
    if idx == 1 || idx == length(f(:))
        disp 'Index at boundaries!'
        warn = 1;
    end
end