%% evaluate performance (i.e. runtime) of different runs

%% setup (user input)
% construct the input file names: prefix, suffix, indexset
pre = 'benchmarks/threads/96_threads_';
suf = '';
ind = [1,2,4,8,12,16];
% what are we looking at
feat = '#threads';

%% run
loadDsets;
require(dsetsRuntime);
nn = length(ind);
times = zeros(nn,1);
for i = 1:nn
    name = [pre num2str(ind(i)) suf];
    readDsets;
    times(i) = runtime_total;
end
plot(ind, times); xlabel(feat); ylabel('runtime [s]'); shg;