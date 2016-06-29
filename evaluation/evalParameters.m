%% list the simulation parameters of a single run

%% setup (user input)
% the filename
name = 'resolution/64_16_5e-3';

%% run
loadDsets;
dsets = dsetsPars;
readDsets;
for i = 1:length(dsetsPars)
    eval(dictNames(dsetsPars{i}))
end