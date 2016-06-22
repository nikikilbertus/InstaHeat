if ~exist('name','var')
    name = 'cutoff/96_16_5e-3';
end
fullname = ['~/Dropbox/Uni/Exercises/12Semester/MAPhysics/data/' name '.h5'];
% fullname = ['~/Data/' name '.h5'];
loadDsets;
if ~exist('dsets', 'var')
    dsets = dsetsAll;
end
require('a', 't', dsets);
% TODO decide whether I want certain indices for the datasets
% create dsetsIndexed in loadDsets to indicate the ones I need to index

for loopvar = 1:length(dsets)
    readDset(dictNames, fullname, dsets{loopvar});
end