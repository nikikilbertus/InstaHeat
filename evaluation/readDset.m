function data = readDset(dictNames, file, dset, I)
    try
        data = h5read(file, ['/' dset]);
    catch
        display(['could not read ' dset])
        return;
    end
    if exist('I','var')
        s = size(data);
        if s(2) > 1
            data = data(:, I);
        else
            data = data(I);
        end
    end
    assignin('caller', dictNames(dset), data);
end