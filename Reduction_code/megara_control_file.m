function [final_data,objectname,star_database]=megara_control_file(objectname,overwrite_fulldata,overwrite_post_reduction,blue_data_chop_value,apply_barycentric_correction,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,order_merge,cross_correlation,teff,logg,vsini,skip_reduction,reduction_folders,ascii_output,fits_output,star_census,FAMIAS_output,rv_measurement,supress_figures,extend_velocities,skip_post_reduction)
warning on
warning('off','backtrace')

final_data=[];
home_dir_megara=which('pointer_file_megara.mat');
home_dir_megara=home_dir_megara(1:end-24);
cd(home_dir_megara)

%setting up directory locations
[raw_data_location,reduced_data_location]=define_directory_locations;

%put single object into a structure if necessary
if isa(objectname,'char') && strcmp(objectname,'musician')==0 && strcmp(objectname,'exocomet')==0
    single_name=objectname;
    clear objectname
    objectname.target1=single_name;
end

%object census
if star_census==1
    found_objectnames=create_observation_summaries(raw_data_location,reduced_data_location,reduction_folders,objectname);
end

object_type=[];
%objectname is predefined target(s)
if isa(objectname,'struct')
    object_type='predefined';
%   %add in wildcard alias names so they are also reduced
%     for l=1:numel(fieldnames(objectname))
%         eval(['single_objectname=objectname.target' num2str(l) ';'])
%         database_obs{l}=single_objectname;
%         if find(ismember(wildc,single_objectname))>0
%             [x,y]=find(ismember(wildc,single_objectname));
%             for n=1:length(wildc(x,:))
%                 eval(['objectname.target' num2str(numel(fieldnames(objectname))+1) '=''' wildc{x(n),2} ''';'])
%                 database_obs{length(database_obs)+1}=wildc{x(n),2};
%                 display(['Alias ' wildc{x(n),2} ' also being reduced. Make sure you skip Processing then check this folder in the reduced data directory for more files.'])
%             end
%         end
%     end   
    for l=1:numel(fieldnames(objectname))
        eval(['single_objectname=objectname.target' num2str(l) ';'])
        display(['Reducing object ' single_objectname])
        try eval(['load(''star_database_' single_objectname '.mat'');'])
            eval(['star_database=star_database_' single_objectname ';'])
        catch
            try eval(['load(''star_database_HD' single_objectname(3:end) '.mat'');'])
                eval(['star_database=star_database_HD' single_objectname(3:end) ';'])
            catch
            star_census=1;
            star_database.exptime=[];
            create_observation_summaries(raw_data_location,reduced_data_location,reduction_folders,objectname)
            end
        end
        [finaldata]=reduce_process_1_target(single_objectname,overwrite_fulldata,overwrite_post_reduction,blue_data_chop_value,apply_barycentric_correction,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,order_merge,cross_correlation,teff,logg,vsini,skip_reduction,rv_measurement,reduction_folders,raw_data_location,reduced_data_location,FAMIAS_output,supress_figures,extend_velocities,skip_post_reduction);
        drawnow
        eval(['final_data.' single_objectname '=finaldata;'])
        if skip_post_reduction==0
            if order_merge==1
                eval(['save(''final_data_' single_objectname ''',''finaldata'')'])
            else
                eval(['save(''final_data_unmerged_' single_objectname ''',''finaldata'')'])
            end
        end
    end
    %objectname is 'musician' or 'tess' or 'exocomet'
elseif strcmp(objectname,'musician') || strcmp(objectname,'tess') || strcmp(objectname,'exocomet')
    object_type='database';
    if strcmp(objectname,'musician')
        load mus_database_cool
        database=mus_database;
    elseif strcmp(objectname,'tess')
        load tess_database
        database=tess_database;
    else
        load exocomet_database
        database=exocomet_database;
    end
    num_objects=0;
    for n=setdiff(1:length(database), [46,48,59])
        single_objectname=database{n,1};
        teff=database{n,2};
        logg=database{n,3};
        vsini=database{n,4};
        num_objects=num_objects+1;
        clear objectname
        eval(['objectname.target' num2str(num_objects) '=''' single_objectname ''';'])
        display(['Reducing object ' single_objectname])
        if star_census==0
            try
                eval(['load(''star_database_' single_objectname '.mat'');'])
                eval(['star_database=star_database_' single_objectname ';'])
            catch
                star_census=1;
                star_database.exptime=[];
            end
        else
            star_database.exptime=[];
        end
        [finaldata]=reduce_process_1_target(single_objectname,overwrite_fulldata,overwrite_post_reduction,blue_data_chop_value,apply_barycentric_correction,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,order_merge,cross_correlation,teff,logg,vsini,skip_reduction,rv_measurement,reduction_folders,raw_data_location,reduced_data_location,FAMIAS_output,supress_figures,extend_velocities,skip_post_reduction);
        drawnow
        eval(['final_data.' single_objectname '=finaldata;'])        
        if skip_post_reduction==0
            if order_merge==1
                eval(['save(''final_data_' single_objectname ''',''finaldata'')'])
            else
                eval(['save(''final_data_unmerged_' single_objectname ''',''finaldata'')'])
            end
        end
    end
    % all stellar frames are reduced
elseif isempty(objectname)
    object_type='all';
    disp(['Reducing all objects'])
    if star_census==0
        found_objectnames=create_observation_summaries(raw_data_location,reduced_data_location,reduction_folders,objectname);
    end
    objectname=found_objectnames;
    fprintf('A total of %.0f objects will be reduced\n',numel(fieldnames(found_objectnames)))
    for l=1:numel(fieldnames(found_objectnames))
        eval(['single_objectname=found_objectnames.target' num2str(l) ';'])
        display(['Reducing object ' single_objectname])
        [finaldata]=reduce_process_1_target(single_objectname,overwrite_fulldata,overwrite_post_reduction,blue_data_chop_value,apply_barycentric_correction,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,order_merge,cross_correlation,teff,logg,vsini,skip_reduction,rv_measurement,reduction_folders,raw_data_location,reduced_data_location,FAMIAS_output,supress_figures,extend_velocities,skip_post_reduction);
        drawnow
        eval(['final_data.' single_objectname '=finaldata;'])
        if skip_post_reduction==0
            if order_merge==1
                eval(['save(''final_data_' single_objectname ''',''finaldata'')'])
            else
                eval(['save(''final_data_unmerged_' single_objectname ''',''finaldata'')'])
            end
        end
    end
    star_database='See individual database files';
end



clearvars -except object_type final_data objectname star_database ascii_output fits_output apply_barycentric_correction rv_measurement cross_correlation
delete('temp*')
cd ../

if ascii_output==1
    write_ascii_output(objectname,final_data,apply_barycentric_correction,rv_measurement,cross_correlation)
end

if fits_output==1
    write_fits_output(objectname,final_data,apply_barycentric_correction,rv_measurement,cross_correlation)
end
