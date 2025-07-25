function[finaldata]=process_1_target(cross_correlation,order_merge,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,apply_barycentric_correction,objectname,teff,logg,vsini,overwrite_fulldata,reduced_data_location,rv_measurement,FAMIAS_output,supress_figures,extend_velocities) 
%PROCESSING MASTER FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(teff) || isempty(logg) || isempty(vsini)
    load mus_database_cool
    flag=0;
    for n=1:length(mus_database)
        dataname=mus_database{n,1};
        if strcmpi(objectname,dataname)
            if isempty(teff)
                teff=mus_database{n,2};
            end
            if isempty(logg)
                logg=mus_database{n,3};
            end
            if isempty(vsini)
                vsini=mus_database{n,4};
            end
            if isempty(mus_database{n,5})
                warning('One or more stellar properties have used average values. Find specific values to improve the CCF')
            end
            flag=1;
        end

    end
    if flag==0
        warning('One or more stellar properties (teff,logg,vsini) have not been set. Find values to input')
    end
end
               
cd([reduced_data_location, '/' ,objectname])

myflag=0;
while(myflag==0)
    [myflag,finaldata]=post_reduction(cross_correlation,order_merge,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,apply_barycentric_correction,objectname,teff,logg,vsini,overwrite_fulldata,rv_measurement,FAMIAS_output,supress_figures,extend_velocities);
end
