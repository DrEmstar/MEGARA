function [finaldata]=reduce_process_1_target(objectname,overwrite_fulldata,overwrite_post_reduction,blue_data_chop_value,apply_barycentric_correction,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,order_merge,cross_correlation,teff,logg,vsini,skip_reduction,rv_measurement,reduction_folders,raw_data_location,reduced_data_location,FAMIAS_output,supress_figures,extend_velocities,skip_post_reduction);

red_flag=0;
%REDUCTION MASTER FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skip_to_post_processing=0;
%picking up any previous reduction files
cd(reduced_data_location)
try movefile('J*_prc.mat','processing_files');
catch
end

%setting up folders in reduced data
files_reduced = dir(reduced_data_location);
dirFlagsr = [files_reduced.isdir];
subFoldersr= files_reduced(dirFlagsr);

%Test for existing reduction of the files in the data directory
flag=0;
for m = 1:length(subFoldersr)
    if strcmpi(subFoldersr(m).name, objectname)
        cd(subFoldersr(m).name)
        try movefile('J*_prc.mat','processing_files');
        catch
        end
        past_reduction=dir('J*.mat');
        flag=1;
        break
    else
        past_reduction=[];
    end
end

if flag==0;
    warning('off','all')
    mkdir(objectname);
    mkdir(objectname,'processing_files')
    warning('on','all')
end

%loading star database
try
    eval(['load(''star_database_' objectname '.mat'');'])
catch
    disp(['No observations found for ' objectname '.'])
    finaldata=[];
    return
end


%Count past reduction files of target found
database_name=whos('star_database_*');
database_name=database_name.name;
eval(['star_database=' database_name ';'])
number_observations=length(star_database.exptime);
past_reduction_list=[];
for u=1:numel(past_reduction)
    past_reduction_list(u)=str2double(past_reduction(u).name(2:end-4));
end
number_observations_reduced=length(past_reduction_list);

if strcmpi(objectname,'dark')
    disp(['Skipping reduction for ' objectname '.'])
    finaldata=[];
    return
end

if strcmpi(objectname,'thorium')
    disp(['Skipping reduction for ' objectname '.'])
    finaldata=[];
    return
end

fprintf('A total of %.0f observations have been found\n',number_observations)
fprintf('A total of %.0f observations have been previously reduced\n',number_observations_reduced)

if skip_reduction==1
    disp(['Skipping reduction for ' objectname '.'])
    
else
    if number_observations==0
        disp(['Nothing to reduce for ' objectname ', please check objectname'])
        finaldata=[];
        return
    end
    %Check number in database
    if number_observations>0
        filename_list=star_database.filename;
        for t=1:numel(filename_list)
            nam_temp=filename_list{t};
            name_list(t)=str2double(nam_temp(2:end-4));
        end
    end
    %Begin reduction if unreduced observations, moves files to target folder
    %with subdirectory with processing files
    
    if numel(reduction_folders)==0 && number_observations>length(past_reduction_list)
        months=star_database.months;
    elseif numel(reduction_folders)~=0
        months=reduction_folders;
    else
        disp('No files to reduce, moving to post-reduction')
        skip_to_post_processing=1;
    end
    if skip_to_post_processing==0;
        for q= 1:length(months)
            cd(raw_data_location)
            cd(months{q})
            fprintf('Reducing Month/Run %.0f of %.0f, This is %s\n',q,length(months),months{q})
            if exist('./file_info.mat', 'file') == 2  % Check only in the current directory
                load('./file_info.mat');  % Load the file from the current directory
            else
                info = display_directory_file_information_nooutput();
                save('file_info','info')
            end
            if numel(info.list)==0
                disp('No files to reduce, moving to next Month/Run')
                finaldata=[];
                return
            end
            start=info.list{1};
            stop=info.list{end};
            startn=str2num(start(2:end));
            stopn=str2num(stop(2:end));
            T=name_list>=startn & name_list<=stopn;
            num_should=length(find(T));
            if isempty(past_reduction)
                num_have=0;
            else
                U=past_reduction_list>=startn & past_reduction_list<=stopn;
                num_have=length(find(U));
            end
            if num_have ~= num_should
                reduce_hercules_frames(objectname,past_reduction,blue_data_chop_value)
                red_flag=1;
                try    eval(['movefile(''J*.mat'',''' reduced_data_location ,'/', objectname ''')'])
                catch
                end
            end
        end
        final_files=dir('J*.mat');
        %barcentric correction all files
        for v=1:numel(final_files)
            load(final_files(v).name)
            nam=final_files(v).name;
            nam(end-3:end)=[];
            eval(['sname=' nam '.sname;'])
            eval(['expdate=' nam '.expdate;'])
            try
                eval(['midtime_tt=' nam '.midtime_tt;'])
            catch
                eval(['midtime_utc=' nam '.midtime_hrsp;'])
                midtime_utc_dt=datetime([expdate(2:end-1) ' ' midtime_utc(2:end-1)],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
                midtime_tt=midtime_utc_dt+seconds(69.184);%%%%%%%%%%%%%%%%
                midtime_tt=['"' datestr(midtime_tt,'HH:MM:SS.FFF') '"'];
                eval([nam '.midtime_utc=midtime_utc;'])
                eval([nam '.midtime_tt=midtime_tt;'])
            end
            eval(['starname=' nam '.starname;'])
            try eval(['bcorr=' nam '.bcorr;'])
            catch
                bcorr=[];
            end
            %collect barycentric correction from C code from HRSP if needed
            if isempty(bcorr) || isnan(bcorr)
                file_name=[nam '.fit'];
                bcorr=str2num(get_barycentric_correction(starname,expdate,midtime_tt,file_name));
                eval([nam '.bcorr=bcorr;'])
                save(nam,nam)
            else
            end
            
            clear(nam)
        end
    end
    cd([reduced_data_location, '/' ,objectname])
    try movefile('J*_prc.mat','processing_files');
    catch
    end
end

cd([reduced_data_location, '/' ,objectname])

try movefile('J*_prc.mat','processing_files');
catch
end

final_files=dir('J*.mat');
if skip_post_reduction==0;
    if overwrite_fulldata==0 && red_flag==0 && size(final_files,1)~=0 && cross_correlation==0 && rv_measurement==0
        eval(['test1=dir(''final_data_' objectname '.mat'');'])
        if numel(test1)>0
            disp('Star previously processed')
            eval(['load(''final_data_' objectname '.mat'');'])
            if supress_figures==0
                figure
                plot(finaldata.fullwave,finaldata.fullint)
                eval(['title(''Full Spectra ' objectname(4:end) ''')'])
            end
        else
            if overwrite_post_reduction==0
                process_data_flag=0;
                try
                    eval(['load(''final_data_' objectname '.mat'');'])
                    if isnan(finaldata.systemic_velocity_direct)
                        process_data_flag=1;
                    end
                catch
                    process_data_flag=1;
                end
                if process_data_flag==1 || overwrite_post_reduction==1
                    finaldata=process_1_target(cross_correlation,order_merge,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,apply_barycentric_correction,objectname,teff,logg,vsini,overwrite_fulldata,reduced_data_location,rv_measurement,FAMIAS_output,supress_figures,extend_velocities)
                end
            end
        end
    elseif size(final_files,1)~=0
        process_data_flag=0;
        try
            eval(['load(''final_data_' objectname '.mat'');'])
        catch
            process_data_flag=1;
        end
        if process_data_flag==1 || overwrite_post_reduction==1 || overwrite_fulldata==1
            finaldata=process_1_target(cross_correlation,order_merge,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,apply_barycentric_correction,objectname,teff,logg,vsini,overwrite_fulldata,reduced_data_location,rv_measurement,FAMIAS_output,supress_figures,extend_velocities)
        end
    else
        finaldata=[];
        process_data_flag=0;
    end
    if  overwrite_post_reduction==0 || process_data_flag==0
        disp('Skipping post-reduction')
    end
else
    try
        eval(['load(''final_data_' objectname '.mat'');'])
    catch
        finaldata=[];
    end
    disp('Skipping processing and post-reduction')
end
