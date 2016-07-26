if ~exist('name','var')
    name = 'cutoff/96_16_5e-3';
end
fullname = ['~/Data/' name '.h5'];
loadDsets;
require('a');

for loopvar = 1:length(dsets)
    if exist('readI','var')
        readDset(dictNames, fullname, dsets{loopvar}, readI);
    else
        readDset(dictNames, fullname, dsets{loopvar});
    end
end