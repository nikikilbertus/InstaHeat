function eq = equalruns(name1, name2)
    name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name1 '.h5'];
    pars1 = getparameters(name);
    name = ['~/Dropbox/Uni/Exercises/11Semester/MAPhysics/data/' name2 '.h5'];
    [pars2, parnames] = getparameters(name);
    eq = true;
    for i = 1:length(pars1)
        if iscell(pars1{i}) && iscell(pars2{i})
            pars1{i} = pars1{i}{:};
            pars2{i} = pars2{i}{:};
        end
        if any(pars1{i} ~= pars2{i})
            eq = false;
            disp(['inequality in component ' parnames{i} ':'])
            disp(pars1{i})
            disp(pars2{i})
        end
    end
end

function [pars, parnames] = getparameters(name)
    parnames = {'/dimension','/time','/gridpoints_internal', ...
        '/gridpoints_output','/mass','/inflaton_mass','/strides_space', ...
        '/strides_time','/tolerances','/max_dt_hubble_fraction','/method', ...
        '/seed','/spatial_bounds_x','/spatial_bounds_y','/spatial_bounds_z'};
    for i = 1:length(parnames)
       pars{i} = h5read(name, parnames{i});
    end
    t = pars{2};
    pars{2} = [t(1) t(end)];
end