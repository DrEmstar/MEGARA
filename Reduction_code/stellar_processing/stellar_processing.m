function stellar_processing(info,allorders,flatdata,jd_ths,fwave,m,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,objectname,past_red)

%%Check for past reductions in raw data folder
% Precompute which .mat files already exist
existingMatFiles = dir('J*.mat');
existingMatNames = {existingMatFiles.name};
existingBaseNames = erase(existingMatNames, '.mat');
if ~isempty(existingBaseNames)
    for n=1:numel(existingBaseNames)
        fprintf('%s.mat exists, skipping this observation\n',existingBaseNames{n})
    end
end

% Only include items that aren't already processed
todoIdx = ~ismember(info.list, existingBaseNames);
filteredList = info.list(todoIdx);
clear existingMatFiles existingBaseNames

disp(['processing stellar images @ ' datestr(now, 'HH:MM:SS.FFF')])

%parallel processing of stellar frames
for s = 1:numel(info.list)%par
    fileobjectname2 = [];
    flag = 0;
    
    [temp, jdexp] = make_temp_variable(info, s);
    
    if strcmpi(temp, 'Stellar')
        stellarimgname = [info.list{s} '.fit'];
        
        % Skip if already reduced
        for ss = 1:length(past_red)
            if strcmpi(past_red(ss).name, [info.list{s} '.mat'])
                fprintf('%s.mat exists, skipping this observation\n', info.list{s});
                flag = 1;
                break;
            end
        end
        if flag > 0
            continue
        end

        % Skip if Julian date is missing
        if isempty(jdexp)
            fprintf('Couldn''t get Julian date for stellar image %s\n', info.list{s});
            continue
        end

        % Find nearest ThAr images
        tmp = sort(jd_ths - jdexp);
        t2 = tmp(tmp < 0);
        indbefore = numel(t2);
        indafter = indbefore + 1;

        if indbefore < 1
            trueindbefore = [];
        else
            trueindbefore = find(round((jd_ths - jdexp) * 1e7) / 1e7 == round(tmp(indbefore) * 1e7) / 1e7);
        end
        if indafter > numel(tmp)
            trueindafter = [];
        else
            trueindafter = find(round((jd_ths - jdexp) * 1e7) / 1e7 == round(tmp(indafter) * 1e7) / 1e7);
        end

        if isempty(trueindbefore) && isempty(trueindafter)
            disp('No ThAr calibration images found')
            continue
        else
            if isempty(trueindbefore)
                trueindbefore = trueindafter;
            end
            if isempty(trueindafter)
                trueindafter = trueindbefore;
            end
        end

        [fwave1, fwave2] = make_fwave_variable(trueindbefore, trueindafter, fwave);

        % Try to get object name
        try
            fileobjectname2 = get_obj_name(stellarimgname);
        catch
            fprintf('Failed to get object name for file %s\n', stellarimgname);
            fileobjectname2 = [];
        end

        if isempty(objectname) || strcmpi(objectname, fileobjectname2)
            fprintf('Reducing stellar spectrum file: %s\n', stellarimgname);
            stellar_frame_processing(stellarimgname, info, allorders, flatdata,fwave1, fwave2,jd_ths(trueindbefore), jd_ths(trueindafter),blue_data_chop_value, flatdata_medfilt, backdata_medfilt,cosmic_plotting, jdexp);
        end
    end

end
