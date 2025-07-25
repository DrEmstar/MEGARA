function [] = reduce_hercules_frames(objectname,past_red,blue_data_chop_value)

%1-program to identify and remove with incomplete/mislabelled data. This runs automatically unless badfiles.mat exisits in the diectory
[change_flag,counter_good]=file_check_routine(objectname);

if change_flag==1
    info = display_directory_file_information_nooutput();
    save('file_info','info')
else
    load file_info.mat
end

fprintf('A total of %.0f stellar frames are suitable for reduction\n',counter_good)

if counter_good==0
    disp('No files to reduce')
    return
end

Bias=[];
plotting=0; % 0 for no, 1 for yes flat field cosmic ray filtering plotting (default=30)
flatdata_medfilt=30; % amount of median smoothing on flat-field data (default=30)
backdata_medfilt=50; % amount of mean smoothing on stellar and flat-field background data (default=50)
cosmic_plotting=0; % whether to plot the cosmic ray filtering process for stellar images (default=0)
numbad_thar=0;%number of thorium files with less than 900 lines found


%2 process the flat-fields by summing all flats, order-tracing and background fitting
[allorders, flatdata, backdata, summed_flat, backfit,save_flag_ff]=flat_processing(info,Bias,blue_data_chop_value,backdata_medfilt,plotting);
if save_flag_ff==1
    eval(['save(''flat_info_blue_' num2str(blue_data_chop_value) ''',''allorders'',''flatdata'',''backdata'',''summed_flat'',''backfit'')'])
end

%3- The following will process all the thar images taken that night
[jd_ths, m, fwave,numbad_thar,save_flag_th] = thar_processing(info, allorders, Bias, flatdata, blue_data_chop_value,numbad_thar);
if save_flag_th==1
    eval(['save(''ThAr_info_blue_' num2str(blue_data_chop_value) ''',''jd_ths'',''m'',''fwave'')'])
end

%4- Process the stellar images one by one using the two thoriums either side of the image (unless there aren't any!).
%include alias names

wildcard_matching=0;
alt_starnames{1}=objectname;
load Mus_known_alias_names
if find(ismember(wildc,objectname))>0
    wildcard_matching=1;
    [x,y]=find(ismember(wildc,objectname));
    for n=1:length(wildc(x,:))
        alt_starnames{n+1}=wildc{x(n),2};
    end
end

for n=1:numel(alt_starnames)
    starname=alt_starnames{n};
    %only process if have files of this name
    process_flag=0;
    for m=1:numel(info.list)
        [type,star] = identify_the_file(info,m);
        if strcmp(star,starname)
            process_flag=1;
        end
    end
    if process_flag==1
    stellar_processing(info,allorders,flatdata,jd_ths,fwave,m,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,starname,past_red);
    end
end



%notify the user that the run has finished
disp(' ')
disp('Run reduced')
if numbad_thar~=0
    warning(' %.0f thoriums may not have reduced correctly!\n',numbad_thar)
end
disp(' ')




