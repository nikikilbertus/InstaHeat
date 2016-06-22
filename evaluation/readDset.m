function data = readDset(dictNames, file, dset, indices)
    try
        data = h5read(file, ['/' dset]);
        if size(data,2) > 1
            data = data';
        end
    catch
        display(['could not read ' dset])
        return;
    end
    if exist('indices','var')
        sizeData = size(data);
        if all((sizeData > 1))
            sizeData
            data = data(indices,:);
        else
            data = data(indices);
        end
    end
    if ~isempty(strfind(dset, 'summary'))
        mkSummary(dictNames, dset, data);
    else
        assignin('caller', dictNames(dset), data);
    end 
end

% TODO: check for constraints, compute rms for phi, psi, ...

function [] = mkSummary(dictNames, dset, data)
    assignin('base', strrep(dictNames(dset),'S','mean'), data(:,1));
    assignin('base', strrep(dictNames(dset),'S','std'), sqrt(data(:,2)));
    assignin('base', strrep(dictNames(dset),'S','min'), data(:,3));
    assignin('base', strrep(dictNames(dset),'S','max'), data(:,4));
    if ~isempty(strfind(dset, 'rho'))
        assignin('base', 'H', sqrt(data(:,1) / 3));
        assignin('base', 'rhorms', sqrt(data(:,2)) ./ data(:,1));
    end
end