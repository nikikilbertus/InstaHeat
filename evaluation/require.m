function [] = require(varargin)
    dsets = joinDsets(varargin{:});
    assignin('caller', 'dsets', dsets);
end