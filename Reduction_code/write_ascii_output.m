function write_ascii_output(objectname,final_data,bcorr_flag,rv_measurement,cross_correlation)

names=fieldnames(objectname);
for m=1:numel(names)
    object=eval(['objectname.' (names{m})]);
    if isempty(eval(['final_data.' object]))==0
        cd(object)
        %Write CCF profiles in dat format
        if cross_correlation==1
            try eval(['load(''ccf_' object '.mat'');'])
                ascii_ccf_file(:,1)=coarse_ccf.velocity;
                for obs=1:numel(coarse_ccf.jd)
                    col=obs+1;
                    ascii_ccf_file(:,col)=coarse_ccf.intensity(obs,:);
                end
                ascii_ccf_filename=['ccf_' object '.dat' ];
                save(ascii_ccf_filename,'ascii_ccf_file','-ascii')
            catch
            end
        end
        %decide if this is order-by-order
        obj_struct=eval(['final_data.' object]);
        obj_struct_names=fieldnames(obj_struct);
        if strncmpi(obj_struct_names{1},'int',3)
            %order-by-order data
            file_list=dir('J*.mat');
            for obs=1:length(file_list)
                ascii_filename=file_list(obs).name(1:end-4);
                num_orders=numel(obj_struct_names)/3;
                for ord=num_orders:-1:1
                    if ord==num_orders
                        %wavelength axis
                        ascii_file(:,1)=eval(['obj_struct.' obj_struct_names{num_orders*3-1}]);
                        %intensity axis
                        ascii_file(:,2)=eval(['obj_struct.' obj_struct_names{num_orders*3-2} '(obs,:)']);
                        %weights axis
                        ascii_file(:,3)=eval(['obj_struct.' obj_struct_names{num_orders*3} '(obs,:)']);
                        %order axis
                        ord_num_real=obj_struct_names{num_orders*3};
                        ord_num_real=str2num(ord_num_real(8:end));
                        ascii_file(:,4)=ord_num_real.*ones(length(eval(['obj_struct.' obj_struct_names{num_orders*3-1}])),1);
                    else
                        %wavelength axis
                        ascii_file_part(:,1)=eval(['obj_struct.' obj_struct_names{ord*3-1}]);
                        %intensity axis
                        ascii_file_part(:,2)=eval(['obj_struct.' obj_struct_names{ord*3-2} '(obs,:)']);
                        %weights axis
                        ascii_file_part(:,3)=eval(['obj_struct.' obj_struct_names{ord*3} '(obs,:)']);
                        %order axis
                        ord_num_real=obj_struct_names{ord*3};
                        ord_num_real=str2num(ord_num_real(8:end));
                        ascii_file_part(:,4)=ord_num_real.*ones(length(eval(['obj_struct.' obj_struct_names{ord*3-1}])),1);
                        ascii_file=[ascii_file ; ascii_file_part];
                        ascii_file_part=[];
                    end
                end
                %summary file
                summary_file{obs,1}=ascii_filename;
                header = get_header([ascii_filename '.fit']);
                headers = get_required_headers_from_header(header);
                header=header';
                load(file_list(obs).name);
                summary_file{obs,2}=eval([ascii_filename '.jd']);
                summary_file{obs,3}=headers.DATE_OBS;
                summary_file{obs,4}=headers.REC_STRT;
                summary_file{obs,5}=headers.EXPTIME;
                summary_file{obs,6}=headers.HERCFWMT;
                summary_file{obs,7}=eval([ascii_filename '.bcorr']);
                summary_file{obs,8}=bcorr_flag;
                summary_file{obs,9}=eval([ascii_filename '.signal_to_noise']);
                ascii_filename=[ascii_filename '_reduced.dat'];
                fid = fopen(ascii_filename,'wt');
                fprintf(fid,'Wavelength Intensity Weight Order \n');
                for l=1:size(ascii_file,1)
                    fprintf(fid,'%.2f %.12f %.12f %.0f\n',ascii_file(l,:));
                end
                fclose(fid);
                %save(ascii_filename,'ascii_file','-ascii')
                ascii_file=[];
                header=[];
                headers=[];
                eval([ascii_filename '=[];'])
            end
            sum_filename=['summary_file_' object '.txt'];
            T = cell2table(summary_file,'VariableNames',{'Filename' 'JD' 'Date' 'Exp_Start' 'Exp_Time' 'FWMT' 'BCorr' 'BCorr_App' 'S/N'});
            writetable(T,sum_filename)
        end
        if strncmpi(obj_struct_names{1},'full',4)
            %full merged spectra
            file_list=dir('J*.mat');
            for obs=1:length(file_list)
                ascii_filename=file_list(obs).name(1:end-4);
                %wavelength axis
                ascii_file(:,1)=obj_struct.fullwave;
                %intensity axis
                ascii_file(:,2)=obj_struct.fullint(obs,:);
                %summary file
                summary_file{obs,1}=ascii_filename;
                header = get_header([file_list(obs).name(1:end-4) '.fit']);
                headers = get_required_headers_from_header(header);
                header=header';
                load(file_list(obs).name);
                summary_file{obs,2}=eval([ascii_filename '.jd']);
                summary_file{obs,3}=headers.DATE_OBS;
                summary_file{obs,4}=headers.REC_STRT;
                summary_file{obs,5}=headers.EXPTIME;
                summary_file{obs,6}=headers.HERCFWMT;
                summary_file{obs,7}=eval([ascii_filename '.bcorr']);
                summary_file{obs,8}=bcorr_flag;
                summary_file{obs,9}=eval([ascii_filename '.signal_to_noise']);
                ascii_filename=[ascii_filename '_reduced.dat'];
                fid = fopen(ascii_filename,'wt');
                fprintf(fid,'Wavelength Intensity \n');
                for l=1:size(ascii_file,1)
                    fprintf(fid,'%.2f %.12f\n',ascii_file(l,:));
                end
                fclose(fid);
                if rv_measurement==1
                    summary_file{obs,10}=obj_struct.systemic_velocity_gauss;
                    try 
                        summary_file{obs,11}=obj_struct.radial_velocity(obs);
                    catch
                        summary_file{obs,11}=NaN;
                    end
                    try
                        summary_file{obs,12}=obj_struct.radial_velocity_error;
                    catch
                        summary_file{obs,12}=NaN;
                    end
                    try
                        summary_file{obs,13}=obj_struct.vsini(obs);
                    catch
                        summary_file{obs,13}=NaN;
                    end
                    summary_file{obs,14}=obj_struct.vsini_error;                 
                else
                end
                ascii_file=[];
                header=[];
                headers=[];
                eval([ascii_filename '=[];'])
            end
            if rv_measurement==0
                sum_filename=['summary_file_' object '.txt'];
                T = cell2table(summary_file,'VariableNames',{'Filename' 'JD' 'Date' 'Exp_Start' 'Exp_Time' 'FWMT' 'BCorr' 'BCorr_App' 'S/N'});
                writetable(T,sum_filename)
            else
                sum_filename=['summary_file_' object '.txt'];
                T = cell2table(summary_file,'VariableNames',{'Filename' 'JD' 'Date' 'Exp_Start' 'Exp_Time' 'FWMT' 'BCorr' 'BCorr_App' 'S/N' 'Syst_Vel' 'RV' 'RV_err' 'Vsini' 'Vsini_err'});
                writetable(T,sum_filename)
            end
        end
        cd ../
    end
end



