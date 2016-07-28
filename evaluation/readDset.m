function data = readDset(dictNames, file, dset, ext, indices)
    try
        data = h5read(file, ['/' dset]);
        if size(data,2) > 1
            data = data';
        end
    catch
        display(['could not read ' dset])
        return;
    end
    loadDsets;
    if exist('indices','var') && contains(dsetsTimeDep,dset)
        sizeData = size(data);
        if all((sizeData > 1))
            data = data(indices,:);
        else
            data = data(indices);
        end
    end
    if ~isempty(strfind(dset, 'summary'))
        mkSummary(dictNames, dset, data, ext);
    else
        assignin('caller', [dictNames(dset) ext], data);
    end 
end

% TODO: check for constraints, compute rms for phi, psi, ...

function [] = mkSummary(dictNames, dset, data, ext)
    assignin('base', strrep(dictNames(dset),'S',['mean' ext]), data(:,1));
    assignin('base', strrep(dictNames(dset),'S',['std' ext]), sqrt(data(:,2)));
    assignin('base', strrep(dictNames(dset),'S',['min' ext]), data(:,3));
    assignin('base', strrep(dictNames(dset),'S',['max' ext]), data(:,4));
    if ~isempty(strfind(dset, 'rho'))
        assignin('base', ['H' ext], sqrt(data(:,1) / 3));
        assignin('base', ['rhorms' ext], sqrt(data(:,2)) ./ data(:,1));
    end
end

function [bool] = contains(dsetsTimeDep, string)
    bool = any(cellfun(@(x) ~isempty(x),strfind(dsetsTimeDep,string)));
end