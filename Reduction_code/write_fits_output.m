function write_fits_output(objectname,final_data,bcorr_flag,rv_measurement,cross_correlation)

import matlab.io.*

names=fieldnames(objectname);
for m=1:numel(names)
    object=eval(['objectname.' (names{m})]);
    if isempty(eval(['final_data.' object]))==0
        cd(object)
        %Write CCF profiles in fits format
        if cross_correlation==1
            fits_ccf_filename=[object '_ccf.fits'];
            try
                fptr = fits.createFile(fits_ccf_filename);
            catch
                delete(fits_ccf_filename)
                fptr = fits.createFile(fits_ccf_filename);
            end
            try eval(['load(''ccf_' object '.mat'');'])
            %image
            ccf_image=[coarse_ccf.velocity' coarse_ccf.intensity'];
            fits.createImg(fptr,'double_img',size(ccf_image));
            fits.writeImg(fptr,ccf_image);
            fits.closeFile(fptr)
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
                fits_filename=file_list(obs).name(1:end-4);
                num_orders=numel(obj_struct_names)/3;
                for ord=num_orders:-1:1
                    fits_image=[];
                    if ord==num_orders
                        %wavelength axis
                        fits_image(:,1)=eval(['obj_struct.' obj_struct_names{num_orders*3-1}]);
                        %intensity axis
                        fits_image(:,2)=eval(['obj_struct.' obj_struct_names{num_orders*3-2} '(obs,:)']);
                        %weights axis
                        fits_image(:,3)=eval(['obj_struct.' obj_struct_names{num_orders*3} '(obs,:)']);
                        %order axis
                        ord_num_real=obj_struct_names{num_orders*3};
                        ord_num_real=str2num(ord_num_real(8:end));
                        fits_image(:,4)=ord_num_real.*ones(length(eval(['obj_struct.' obj_struct_names{num_orders*3-1}])),1);
                    else
                        %wavelength axis
                        fits_file_part(:,1)=eval(['obj_struct.' obj_struct_names{ord*3-1}]);
                        %intensity axis
                        fits_file_part(:,2)=eval(['obj_struct.' obj_struct_names{ord*3-2} '(obs,:)']);
                        %weights axis
                        fits_file_part(:,3)=eval(['obj_struct.' obj_struct_names{ord*3} '(obs,:)']);
                        %order axis
                        ord_num_real=obj_struct_names{ord*3};
                        ord_num_real=str2num(ord_num_real(8:end));
                        fits_file_part(:,4)=ord_num_real.*ones(length(eval(['obj_struct.' obj_struct_names{ord*3-1}])),1);
                        fits_image=[fits_image ; fits_file_part];
                        fits_file_part=[];
                    end
                end
                %header
                header_file{obs,1}=fits_filename;
                header = get_header([file_list(obs).name(1:end-4) '.fit']);
                headers = get_required_headers_from_header(header);
                header=header';
                load(file_list(obs).name);
                header_file{obs,2}=eval([fits_filename '.jd']);
                header_file{obs,3}=headers.DATE_OBS;
                header_file{obs,4}=headers.REC_STRT;
                header_file{obs,5}=headers.EXPTIME;
                header_file{obs,6}=headers.HERCFWMT;
                header_file{obs,7}=eval([fits_filename '.bcorr']);
                header_file{obs,8}=bcorr_flag;
                header_file{obs,9}=eval([fits_filename '.signal_to_noise']);
                fits_filename=[fits_filename '_reduced.fits'];
                %write to fits file
                try
                    fptr = fits.createFile(fits_filename);
                catch
                    delete(fits_filename)
                    fptr = fits.createFile(fits_filename);
                end
                %image
                fits.createImg(fptr,'double_img',size(fits_image));
                fits.writeImg(fptr,fits_image);
                %header
                fits.writeKey(fptr,'Filename',header_file{obs,1},'Original filename');
                fits.writeKey(fptr,'JD',header_file{obs,2},'Reduced, heliocentric, mid-obs JD');
                fits.writeKey(fptr,'Date',header_file{obs,3},'Observation Date');
                fits.writeKey(fptr,'Exp_Time',header_file{obs,4},'Exposure duration');
                fits.writeKey(fptr,'Start',header_file{obs,5},'Exposure start time UTC');
                fits.writeKey(fptr,'FWMT',header_file{obs,6},'Flux-weighted exposure midtime UTC');
                fits.writeKey(fptr,'BCorr',header_file{obs,7},'Barycentric correction (velocity)');
                fits.writeKey(fptr,'Bcorr_App',header_file{obs,8},'Barycentric correction applied?');
                fits.writeKey(fptr,'S/N',header_file{obs,9},'Signal-to-noise');
                fits.closeFile(fptr)
                fits_file=[];
                header=[];
                headers=[];
                eval([fits_filename '=[];'])
            end
        end
        if strncmpi(obj_struct_names{1},'full',4)
            %full merged spectra
            file_list=dir('J*.mat');
            for obs=1:length(file_list)
                fits_filename=file_list(obs).name(1:end-4);
                %wavelength axis
                fits_image(:,1)=obj_struct.fullwave;
                %intensity axis
                fits_image(:,2)=obj_struct.fullint(obs,:);
                %header
                header_file{obs,1}=fits_filename;
                header = get_header([file_list(obs).name(1:end-4) '.fit']);
                headers = get_required_headers_from_header(header);
                header=header';
                load(file_list(obs).name);
                header_file{obs,2}=eval([fits_filename '.jd']);
                header_file{obs,3}=headers.DATE_OBS;
                header_file{obs,4}=headers.REC_STRT;
                header_file{obs,5}=headers.EXPTIME;
                header_file{obs,6}=headers.HERCFWMT;
                header_file{obs,7}=eval([fits_filename '.bcorr']);
                header_file{obs,8}=bcorr_flag;
                header_file{obs,9}=eval([fits_filename '.signal_to_noise']);
                if rv_measurement==1
                    header_file{obs,10}=obj_struct.systemic_velocity;
                    header_file{obs,11}=obj_struct.radial_velocity(obs);
                    header_file{obs,12}=obj_struct.radial_velocity_error;
                    header_file{obs,13}=obj_struct.vsini(obs);
                    header_file{obs,14}=obj_struct.vsini_error;
                end
                fits_filename=[fits_filename '_reduced.fits'];
                %write to fits file
                try
                    fptr = fits.createFile(fits_filename);
                catch
                    delete(fits_filename)
                    fptr = fits.createFile(fits_filename);
                end
                %image
                fits.createImg(fptr,'double_img',size(fits_image));
                fits.writeImg(fptr,fits_image);
                %header
                fits.writeKey(fptr,'Filename',header_file{obs,1},'Original filename');
                fits.writeKey(fptr,'JD',header_file{obs,2},'Reduced, heliocentric, mid-obs JD');
                fits.writeKey(fptr,'Date',header_file{obs,3},'Observation Date');
                fits.writeKey(fptr,'Exp_Time',header_file{obs,4},'Exposure duration');
                fits.writeKey(fptr,'Start',header_file{obs,5},'Exposure start time UTC');
                fits.writeKey(fptr,'FWMT',header_file{obs,6},'Flux-weighted exposure midtime UTC');
                fits.writeKey(fptr,'BCorr',header_file{obs,7},'Barycentric correction (velocity)');
                fits.writeKey(fptr,'Bcorr_App',header_file{obs,8},'Barycentric correction applied?');
                fits.writeKey(fptr,'S/N',header_file{obs,9},'Signal-to-noise');
                if rv_measurement==1
                    fits.writeKey(fptr,'Syst_Vel',header_file{obs,9},'Systemic Velocity');
                    fits.writeKey(fptr,'RV',header_file{obs,10},'Radial Velocity');
                    fits.writeKey(fptr,'RV_err',header_file{obs,11},'Radial Velocity error');
                    fits.writeKey(fptr,'Vsini',header_file{obs,12},'Projected Rotational Velocity vsini');
                    fits.writeKey(fptr,'Vsini_err',header_file{obs,13},'Vsini error');
                end
                fits.closeFile(fptr)
                fits_file=[];
                header=[];
                headers=[];
                eval([fits_filename '=[];'])
            end
        end
        cd ../
    end
end
