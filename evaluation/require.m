function [] = require(varargin)
    if evalin('caller','exist(''dsets'',''var'')')
        dsets = joinDsets(varargin{:}, evalin('caller','dsets'));
    else
        dsets = joinDsets(varargin{:});
    end
    assignin('caller', 'dsets', dsets);
end