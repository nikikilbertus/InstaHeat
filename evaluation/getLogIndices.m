function [I] = getLogIndices(N, nout)
    nouttmp = 0; i = 0;
    while (nouttmp < nout) 
        I = unique(round(logspace(0.1,log10(N),nout+i)));
        nouttmp = length(I);
        i = i+1;
    end
end