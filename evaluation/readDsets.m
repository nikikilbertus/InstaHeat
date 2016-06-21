name = 'cutoff/96_16_5e-3_fil';
name = ['~/Dropbox/Uni/Exercises/12Semester/MAPhysics/data/' name '.h5'];
loadDsets;
if ~exist('dsets', 'var')
    dsets = dsetsAll;
end

% TODO decide whether I want certain indices for the datasets
% create dsetsIndexed in loadDsets to indicate the ones I need to index

for i = 1:length(dsets)
    readDset(dictNames, name, dsets{i});
end