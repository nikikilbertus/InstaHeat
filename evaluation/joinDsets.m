function [dset] = joinDsets(varargin)
    narginchk(1,inf);
    n = nargin;
    dset = cellify(varargin{n});
    n = n - 1;
    while n
        dset = [dset, cellify(varargin{n})];
        n = n - 1;
    end
    dset = unique(dset);
end

function [cell] = cellify(in)
    if ~iscell(in)
        in = namify(in);
        cell = {in};
    else
    	cell = in;
    end
end

function [cell] = namify(in)
    loadDsets;
    res = find(strcmp(in, dsetsAll), 1);
    if isempty(res)
        res = find(strcmp(in, dsetsNames), 1);
        if isempty(res)
            error(['dataset ' in ' does not exist']);
        end
        cell = dictNamesInv(in);
    else
       cell = in; 
    end
end