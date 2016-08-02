%% evaluate performance (i.e. runtime) of different runs

%% setup (user input)
% construct the input file names: prefix, suffix, indexset
pre = 'benchmarks/splitthreads/max32/';
suf = '';
ind1 = [1,4,8,12,16,20,24,28,32];
ind2 = ind1;
% what are we looking at
feat1 = '#threads_omp';
feat2 = '#threads_fftw';

%% run
loadDsets;
require(dsetsRuntime);
nn1 = length(ind1);
nn2 = length(ind2);
times = zeros(nn1,nn2);
for i = 1:nn1
    for j = 1:nn2
        name = [pre num2str(ind1(i)) '_' num2str(ind2(j)) suf];
        readDsets;
        times(i,j) = runtime_total;
    end
end
% plot(ind, times); xlabel(feat); ylabel('runtime [s]'); shg;
% surf(ind1,ind2,times); xlabel(feat1); ylabel(feat2); zlabel('runtime [s]'); shg;
bar3(times); xlabel(feat1); ylabel(feat2); zlabel('runtime [s]'); shg;