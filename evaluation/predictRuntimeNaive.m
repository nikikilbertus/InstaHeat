%% roughly predict how much longer a run will take on a different grid
% Nold and Nnew are the total number of gridpoints (e.g. N^3 in 3 dimensions)
function [factor, runtime] = predictRuntimeNaive(Nold, Nnew, runtimeOld)
    factor = Nnew^3*log(Nnew^3) / (Nold^3*log(Nold^3));
    if exist('runtimeOld','var')
        runtime = factor * runtimeOld;
    else
        runtime = -1;
    end
end