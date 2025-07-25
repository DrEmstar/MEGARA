function [passback_objectname]=create_observation_summaries(raw_data_location,reduced_data_location,reduction_folders,reduction_objectname)

build_database=0;%If this is turned on the database is built from all stellar files in the directory
wildcard_matching=0;
database_build_counter=0;
load Mus_known_alias_names
if strcmp(reduction_objectname,'musician') || strcmp(reduction_objectname,'tess')
    if strcmp(reduction_objectname,'musician')
        load mus_database_cool.mat
        database_obs=mus_database;
    end
    if strcmp(reduction_objectname,'tess')
        load tess_database.mat
        database_obs=tess_database;
    end
    wildcard_matching=1;
    count=zeros(length(database_obs),1);
    count2=zeros(length(database_obs),1);
    port_4_count=zeros(length(database_obs),1);
    fibre_3_count=zeros(length(database_obs),1);
end

if isa(reduction_objectname,'struct')
    for l=1:numel(fieldnames(reduction_objectname)) 
        eval(['single_objectname=reduction_objectname.target' num2str(l) ';'])
        database_obs{l}=single_objectname;
        count=zeros(length(database_obs),1);
        count2=zeros(length(database_obs),1);
        port_4_count=zeros(length(database_obs),1);
        fibre_3_count=zeros(length(database_obs),1);
    end
end

if isempty(reduction_objectname)
    database_obs=[];
    build_database=1;
    count=zeros(1000,1);
    count2=zeros(1000,1);
    port_4_count=zeros(1000,1);
    fibre_3_count=zeros(1000,1);
    database_build_counter=0;
end

% Get a list of all observing folders
files = dir(raw_data_location);
dirFlags = [files.isdir];
subFolders = files(dirFlags);

cd(raw_data_location)

%set up reduction folder subsets
if numel(reduction_folders)>0
    subFolders_selected.name=[];
    c=2;
    cell_subFolders=struct2cell(subFolders);
    cell_subFolders=cell_subFolders(1,:);
    for n=1:numel(reduction_folders)
        tf=strcmp(reduction_folders(n),cell_subFolders);
        hits=find(tf);
        if numel(hits)==0
            disp(['Folder ' reduction_folders{n} ' does not exist in the raw data directory, please check name and directory location'])
        else
            c=c+1;
            subFolders_selected(c).name=subFolders(hits).name;
        end
    end
    subFolders=subFolders_selected;
end

% Collect all the file_information in all folders
port_flag=0;
fibre_flag=0;

for n= 3:length(subFolders)
    fprintf('Running census on Month/Run %.0f of %.0f\n',n-2,length(subFolders)-2)
    cd(raw_data_location)
    cd(subFolders(n).name)
    flag2=0;
    fibre_error_num=0;
    try load file_info
    catch
        info = display_directory_file_information_nooutput();
        save('file_info','info')
    end
    
    %skip checking again if reduction already happened here
    if isfolder('removed_files') || isfolder('Removed files')
    else
        %general bad file check
        badfile_check_general;
        file_list=dir('J*.fit');
        
        for t=1:numel(file_list)
            eval(['header = get_header(''' file_list(t).name ''');'])
            temp_head = get_required_headers_from_header(header);
            port_value = temp_head.PORT;
            fibre_value = temp_head.HERCFIB;
            if port_value > 1
                port_flag=port_flag+1;
                if port_flag==1
                    mkdir 'four_port_files'
                end
                movefile(file_list(t).name,'four_port_files')
            end
            if fibre_value==3
                fibre_flag=fibre_flag+1;
                if fibre_flag==1
                    mkdir 'fibre_3_files'
                end
                movefile(file_list(t).name,'fibre_3_files')
            elseif fibre_value==1
                continue
                %leave fibre error files in dir
            else
                fibre_error_num=fibre_error_num+1;
                eval(['fibre_error_list.name' num2str(fibre_error_num) '=''' file_list(t).name ''';'])
            end
        end
        info = display_directory_file_information_nooutput();
        save('file_info','info')
    end
    
    for s=1:length(info.list)
        [type,star] = identify_the_file(info,s); %I.D. the file
        if build_database==1
            if strcmpi(type, 'Stellar') && sum(strcmpi(database_obs, star))==0
                database_build_counter=database_build_counter+1;
                database_obs{database_build_counter}=star;
                
                eval(['reduction_objectname.target' num2str(database_build_counter) '=''' star ''';'])
            end
        end
        for u=1:length(database_obs)
            objectname=database_obs{u};
            if contains(objectname,'-')
                location=strfind(objectname,'-');
                objectname(location)='_';
            end
            flag2(u,s)=0;
            objectname2=['HD' objectname(4:end)];
            if strcmpi(type, 'Stellar') && strcmpi(star,objectname)
                count(u)=count(u)+1;
                eval(['star_database_' objectname '.filename{count(u)}=''' info.list{s} '.fit'';'])
                eval(['star_database_' objectname '.date{count(u)}=eval([''info.'' info.list{s} ''.DATE_OBS'']);'])
                eval(['star_database_' objectname '.jd{count(u)}=eval([''info.'' info.list{s} ''.MJD_OBS'']);'])
                eval(['star_database_' objectname '.exptime(count(u))=eval([''info.'' info.list{s} ''.EXPTIME'']);'])
                flag2(u,s)=flag2(u,s)+1;
                if sum(flag2(u,:))==1
                    count2(u)=count2(u)+1;
                    eval(['star_database_' objectname '.months{count2(u)}=subFolders(n).name;'])
                end
            elseif strcmpi(type, 'Stellar') && strcmpi(star,objectname2)
                count(u)=count(u)+1;
                eval(['star_database_' objectname '.filename{count(u)}=''' info.list{s} '.fit'';'])
                eval(['star_database_' objectname '.date{count(u)}=eval([''info.'' info.list{s} ''.DATE_OBS'']);'])
                eval(['star_database_' objectname '.jd{count(u)}=eval([''info.'' info.list{s} ''.MJD_OBS'']);'])
                eval(['star_database_' objectname '.exptime(count(u))=eval([''info.'' info.list{s} ''.EXPTIME'']);'])
                flag2(u,s)=flag2(u,s)+1;
                if sum(flag2(u,:))==1
                    count2(u)=count2(u)+1;
                    eval(['star_database_' objectname '.months{count2(u)}=subFolders(n).name;'])
                end
            end
        end
    end
    %wildcard matching for musician targets
    if wildcard_matching==1
        for w=1:length(wildc)
            wildn=wildc{w,1};
            countw=0;
            for n=1:numel(wildc(w,:))-1
                wilda=wildc{w,n+1};
                if exist(['star_database_' wilda])==1
                    if exist(['star_database_' wildn])==1
                        countw=length(eval(['star_database_' wildn '.filename']));
                    end
                    for x=1:length(eval(['star_database_' wilda '.filename']))
                        countw=countw+1;
                        eval(['star_database_' wildn '.filename{countw}=star_database_' wilda '.filename{x};']);
                        eval(['star_database_' wildn '.date{countw}=star_database_' wilda '.date{x};']);
                        eval(['star_database_' wildn '.jd{countw}=star_database_' wilda '.jd{x};']);
                        eval(['star_database_' wildn '.exptime(countw)=star_database_' wilda '.exptime(x);']);
                    end
                    for a=1:numel(database_obs)
                        namst=database_obs{a};
                        if strcmpi(wildn,namst)==1
                            ualt=a;
                            count(ualt)=count(ualt)+x;
                            count2(ualt)=count2(ualt)+1;
                            eval(['star_database_' wildn '.months{count2(ualt)}=star_database_' wilda '.months{end};'])
                        end
                        if strcmpi(wilda,namst)==1
                            ualt=a;
                            count(ualt)=0;
                        end
                    end
                end
                clear(['star_database_' wilda])
            end
        end
    end
end

% if port_flag>0
%     fprintf('A total of %.0f observations are 4-port readout and not reduced\n',port_flag)
% end
% if fibre_flag>0
%     fprintf('A total of %.0f observations are fibre 3 and not reduced\n',fibre_flag)
% end
% if fibre_error > 0;
%     fprintf(' A total of %.0f observations have no recorded fibre position and may not reduce\n',fibre_error)
% end

cd(reduced_data_location)
try cd('Observation Summaries')
catch
    mkdir('Observation Summaries')
    addpath('Observation Summaries')
    cd('Observation Summaries')
end
for v=1:length(database_obs)
    master_file{v,1}=database_obs{v};
    ex=exist(['star_database_' database_obs{v}]);
    if ex==1
        master_file{v,2}=length(eval((['star_database_' database_obs{v} '.filename'])));
        master_file{v,3}=port_4_count(v);
        master_file{v,4}=fibre_3_count(v);
        %clean up months
        eval(['star_database_' database_obs{v} '.months=unique(star_database_' database_obs{v} '.months);'])
        save(['star_database_' database_obs{v}],['star_database_' database_obs{v}])
        
    end    
end

save('observation_summary_count','master_file')

passback_objectname=reduction_objectname;

