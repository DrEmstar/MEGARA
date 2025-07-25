function [myflag,finaldata] = post_reduction(cross_correlation,order_merge,apply_continuum_fit,manual_continuum_fit,auto_continuum_fit,apply_barycentric_correction,objectname,teff,logg,vsini,overwrite_fulldata,rv_measurement,FAMIAS_output,supress_figures,extend_velocities)

%DEFINITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigm=0.5;wstart=3800;wstop=8000;
%systemic velocity
shift=0;
remind=[];
firstorder=76;
lastorder=146;
step=0.02;

radial_velocity2=NaN;
vsini_gauss2=NaN;

weak_limit1=1; %lines < XX mA in strength removed from ccf (default 1)
weak_limit2=5; %lines < YY mA in strength removed from lsd delta function ccf (default 5)

%list files
list=(dir('J*.mat'));
% for m=length(remind):-1:1
%     list(remind(m),:)=[];
% end
cd processing_files/
listp=(dir('*.mat'));
% for m=length(remind):-1:1
%     listp(remind(m),:)=[];
% end
cd ..

test_frame=load(list(1).name);
eval(['test_data=test_frame.' list(1).name(1:end-4) ';'])

%checks maximum order that has been reduced
g=lastorder;
order_exist_test=0;
while order_exist_test==0
    eval(['order_exist_test=isfield(test_data,''order_' num2str(g) ''');'])
    g=g-1;
end
lastorder=g+1;

test1=dir('final_data*.mat');

removed_files_counter=0;
if numel(test1)==0 || (numel(test1)>0 && overwrite_fulldata==1)
    %STAGE I
    fprintf('Merging orders of all frames\n')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Finding all the max/mins for orders
    maxwave=ones(size(firstorder:lastorder))*1e6;
    minwave=zeros(size(firstorder:lastorder));
    for obsnum=1:length(list)
        load(list(obsnum).name)
        nam=list(obsnum).name;
        nam(end-3:end)=[];
        cntt=0;
        for ord=firstorder:lastorder
            cntt=cntt+1;
            eval(['wave=' nam '.order_' num2str(ord) '(:,1);'])
            if min(wave)>minwave(cntt)
                minwave(cntt)=min(wave);
            end
            if max(wave)<maxwave(cntt)
                maxwave(cntt)=max(wave);
            end
        end
        clear(nam)
    end
    
    %Creating the order matrices in the structure array
    for obsnum=1:length(list)
        eval(['load ' list(obsnum,1).name]);
        cd processing_files/
        eval(['load ' listp(obsnum,1).name]);
        cd ..
        file=list(obsnum,1).name(1:8);
        cnt=0;
        for ord=firstorder:lastorder
            cnt=cnt+1;
            if obsnum==1
                eval(['data.wave' num2str(ord) '=rounding(minwave(cnt)+2,0.02):0.02:rounding(maxwave(cnt)-2,0.02);'])
            end
            if apply_barycentric_correction==1;
                try eval(['bcorr=' list(obsnum).name(1:8) '.bcorr;'])
                    bcorr=bcorr*-1;
                catch
                    bcorr=[];
                end
                %calculate barycentric correction if it is missing
                if isempty(bcorr)==1
                    file_name=[file '.fit'];
                    eval(['expdate=' file '.expdate;'])
                    try eval(['midtime_tt=' file '.midtime_tt;'])
                    catch
                        eval(['midtime_utc=' file '.midtime_hrsp;'])
                        midtime_utc_dt=datetime([expdate(2:end-1) ' ' midtime_utc(2:end-1)],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
                        midtime_tt=midtime_utc_dt+seconds(69.184);%%%%%%%%%%%%%%%%update as needed
                        midtime_tt=['"' datestr(midtime_tt,'HH:MM:SS.FFF') '"'];
                        eval([file '.midtime_utc=midtime_utc;'])
                        eval([file '.midtime_tt=midtime_tt;'])
                    end
                    bcorr=str2num(get_barycentric_correction(objectname,expdate,midtime_tt,file_name));
                    eval([file '.bcorr=bcorr;'])
                    save(file,file)
                    bcorr=bcorr*-1;
                end
                if isempty(bcorr)==1
                    eval(['tmpwave=' file '.order_' num2str(ord) '(:,1);'])
                    bcorr=0;
                else
                    eval(['tmpwave=shift_velocity(' file '.order_' num2str(ord) '(:,1),' num2str(bcorr) ');'])
                end
                eval(['data.bcorr(' num2str(obsnum) ')=bcorr;'])
            else
                eval(['tmpwave=' file '.order_' num2str(ord) '(:,1);'])
            end
            eval(['int=' file '.order_' num2str(ord) '(:,2);'])
            eval(['iind1=' file '_prc.order_cutpixels.order_' num2str(ord) '.ind1;'])
            eval(['iind2=' file '_prc.order_cutpixels.order_' num2str(ord) '.ind2;'])
            eval(['weights=' file '_prc.order_nff_' num2str(ord) '(iind1:iind2);'])
            eval(['data.int' num2str(ord) '(obsnum,:)=spline(tmpwave,int,data.wave' num2str(ord) ');'])
            eval(['data.weights' num2str(ord) '(obsnum,:)=spline(tmpwave,weights,data.wave' num2str(ord) ');'])
            eval(['jd(obsnum)=' list(obsnum).name(1:8) '.jd;'])
            try
                eval(['signal_to_noise(obsnum)=' list(obsnum).name(1:8) '.signal_to_noise;'])
            catch
                signal_to_noise(obsnum)=NaN;
            end
        end
        if exist('bcorr','var')==1;
            allbcorr(obsnum)=bcorr;
        end
        eval(['clear ' file])
        eval(['clear ' file '_prc'])
    end
    data.jd=jd;
    for ind3=lastorder:-1:firstorder
        eval(['wave=data.wave' num2str(ind3) ';'])
        %Checking for empty fields
        if isempty(wave)
            eval(['data=rmfield(data,''wave' num2str(ind3) ''');'])
            eval(['data=rmfield(data,''int' num2str(ind3) ''');'])
            eval(['data=rmfield(data,''weights' num2str(ind3) ''');'])
            lastorder=ind3-1;
        end
    end
    
    
    if exist('allbcorr','var')==1;
        data.bcorr=allbcorr;
    end

    
    %STAGE II
    fprintf('Normalising all orders\n')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This section removes order to order differences and fits a polynomial to the data for a coarse normalisation.
    for ind3=firstorder:lastorder
        small_flag=0;
        %fprintf('processing order %.0f of %.0f\n',ind3-firstorder+1,lastorder-firstorder+1)
        eval(['int=data.int' num2str(ind3) ';'])
        eval(['wave=data.wave' num2str(ind3) ';'])
        if size(int,2)>1000
            int=int(:,50:end-50);
            wave=wave(50:end-50);
            small_flag=1;
        end
        int(int>2)=2;
        intx=medfilt1(int',71)';
        intxx=intx(:,71:end-70);
        wavex=wave(71:end-70);
        medianspec=median(intxx);
        for ind1=1:numel(int(:,1))
            temp=intxx(ind1,:);
            temp_m=temp./medianspec;
            warning off
            [p,S,mu]=polyfit(wavex,temp_m,5);
            warning on
            thefit=polyval(p,wave,S,mu);
            into=int(ind1,:)./thefit;
            eval(['data2.int' num2str(ind3) '(ind1,:)=into;'])
        end
        eval(['data2.wave' num2str(ind3) '=wave;'])
        if small_flag==1
            eval(['data2.weights' num2str(ind3) '=data.weights' num2str(ind3) '(:,50:end-50);'])
        else
            eval(['data2.weights' num2str(ind3) '=data.weights' num2str(ind3) ';'])
        end
        eval(['data.int' num2str(ind3) '=[];'])
        eval(['data.wave' num2str(ind3) '=[];'])
        eval(['data.weights' num2str(ind3) '=[];'])
    end
    clear data
    %save('data2','data2')
    %OPTION: CFIT ROUTINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if apply_continuum_fit==1
        fprintf('CONTINUUM FITTING\n')
        if manual_continuum_fit==0 && auto_continuum_fit==0
            customdatacf=eval(['''datacf_' objectname '''']);
            try load(customdatacf)
                fprintf('Loaded previous continuum fit\n')
            catch
                fprintf('No previous continuum fit. Opening program for manual fit\n')
                fprintf('Try to keep the function smooth and fit the top sides of the line profiles.\n')
                datacf=manual_continuum_fitting(teff,logg,vsini,sigm,wstart,wstop,shift,firstorder,lastorder,data2,objectname);
            end
        end
        if manual_continuum_fit==1 
            fprintf('Create new manual continuum fit. Opening program for manual fit\n')
            fprintf('Try to keep the function smooth and fit the top sides of the line profiles.\n')
            datacf=manual_continuum_fitting(teff,logg,vsini,sigm,wstart,wstop,shift,firstorder,lastorder,data2,objectname);
            manual_continuum_fit=0;
        end
        if auto_continuum_fit==1
            load datacf_template
            if lastorder>150
                lastorder=150;
            end
            fprintf('Loaded template continuum fit\n')
        end
        %apply continuum fit
        warning off
        for ord=[firstorder:lastorder]
            eval(['int=data2.int' num2str(ord) ';'])
            eval(['wave=data2.wave' num2str(ord) ';'])
            eval(['fiti=datacf.int' num2str(ord) ';'])
            eval(['fitw=datacf.w' num2str(ord) ';'])
            if wave(1) ~= fitw(1) || numel(wave) ~= numel(fitw)
                fitsp=spline(fitw,fiti,wave);
                avgspec=mean(fiti);
                S=abs(avgspec-fitsp)>0.1;
                Slow=logical([S(1:end/2) zeros(1,round(numel(S)/2))]);
                Shigh=logical([zeros(1,round(numel(S)/2)) S((end/2+1):end)]);
                fitsp(Slow)=fiti(1);
                fitsp(Shigh)=fiti(length(fiti));
                fitsp2 = smooth(fitsp,100);
                eval(['data2.int' num2str(ord) '=data2.int' num2str(ord) './repmat(fitsp2'',size(data2.int' num2str(ord) ',1),1);'])
            else
                eval(['data2.int' num2str(ord) '=data2.int' num2str(ord) './repmat(fiti,size(data2.int' num2str(ord) ',1),1);'])
            end
        end
        warning on
        
        %1D cosmic ray removal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %final cosmic ray removal
        
        fprintf('Removing Cosmic Rays\n')
        for ord=[firstorder:lastorder]
            %fprintf('processing order %.0f\n',ord)
            eval(['int=data2.int' num2str(ord) ';'])
            temp2=medfilt2(int-1,[7 31])+1;
            temp=int-repmat(median(int),[size(temp2,1) 1]);
            S=abs(temp)>6*repmat(std(temp,[],2),1,size(temp,2));
            int(S)=temp2(S);
            int=int(:,31:end-31);
            eval(['data2.int' num2str(ord) '=int;'])
            eval(['data2.wave' num2str(ord) '=data2.wave' num2str(ord) '(31:end-31);'])
            eval(['data2.weights' num2str(ord) '=data2.weights' num2str(ord) '(:,31:end-31);'])
            clear temp*
        end
        
        %OPTION:Merge orders
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Checking for empty fields
        for ind3=lastorder:-1:firstorder
            eval(['wave=data2.wave' num2str(ind3) ';'])
            if isempty(wave)
                eval(['data2=rmfield(data2,''wave' num2str(ind3) ''');'])
                eval(['data2=rmfield(data2,''int' num2str(ind3) ''');'])
                eval(['data2=rmfield(data2,''weights' num2str(ind3) ''');'])
                lastorder=ind3-1;
            end
        end
        
        if order_merge==1;
            fprintf('Merging orders\n')
            %merging orders
            [fullwave,fullint] = combining_orders_2Dint(data2,step,firstorder,lastorder);
            if isempty(remind)==1
                medspec=median(fullint);
                standdev=std(fullint);
                count=0;
                for n=1:size(fullint,1)
                    testspec=fullint(n,:);
                    SS=testspec>5*standdev+medspec;
                    if length(find(SS))>500
                        count=count+1;
                        remind(count)=n;
                    end
                end
                if isempty(remind)==0
                    total_removed=numel(remind)+removed_files_counter;
                    for nn=1:numel(remind)
                        %make removed file directories
                        warning off
                        if nn==1
                            try mkdir 'removed_reduced_files'
                            catch
                            end
                            try mkdir 'processing_files/removed_reduced_files'
                            catch
                            end
                            warning on
                        end
                        %move removed files
                        try movefile(list(remind(nn)).name,'removed_reduced_files');
                        catch
                        end
                        prc_name=[list(remind(nn)).name(1:8) '_prc.mat'];
                        cd processing_files/
                        try movefile(prc_name,'removed_reduced_files');
                        catch
                        end
                        cd ../
                    end
                    fprintf('A total of %.0f observations have been removed\n',total_removed)
                    fprintf('Reprocessing improved dataset\n')
                    finaldata=[];
                    myflag=0;
                    return
                else
                    total_removed=removed_files_counter;
                    fprintf('A total of %.0f observations have been removed\n',total_removed)
                end
            end
            finaldata.fullwave=fullwave';
            finaldata.fullint=fullint;
            finaldata.jd=jd;
            eval(['save(''final_data_' objectname '.mat'',''finaldata'')'])
        else
            finaldata=data2;
            clear data2;
            eval(['save(''final_data_' objectname '.mat'',''finaldata'')'])
        end
    end
else
    disp('previous final_data loaded')
    load(test1.name)
    fullwave=finaldata.fullwave;
    fullint=finaldata.fullint;
    jd=finaldata.jd;
end

%OPTION:Cross-correlation Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cross_correlation==1 || rv_measurement==1
    fprintf('Cross-correlating spectra\n')
    %Find Nans and set to 1
    G=isnan(fullint);
    fullint(G)=1;
    
    [synthwave, synthint, wavelengths, ele, elem, widths] = make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop);
    w=wavelengths;
    wi=widths;
    %cut out Hydrogen and telluric regions from delta function list. The below are the regions to keep.
    %     %general list
    S1=w>4260 & w<4320;
    S2=w>4380 & w<4830;
    S3=w>4895 & w<5865;
    S4=w>5982 & w<6274;
    S5=w>6326 & w<6465;
    S6=w>6600 & w<6858;
    S7=w>7415 & w<7500;
    
    W=cat(1,w(S1),w(S2),w(S3),w(S4),w(S5),w(S6),w(S7));
    Wi=cat(1,wi(S1),wi(S2),wi(S3),wi(S4),wi(S5),wi(S6),wi(S7));
    
    %remove very weak lines
    SS=Wi<weak_limit1;
    W(SS)=[];
    Wi(SS)=[];
    
    %do a ccf with the standard synthetic spectrum wavelengths and widths
    [v,vel,LSD,CCF,dfn_auto,wcw] = do_ccf_from_spec_wavelengths_and_widths(fullwave,fullint,W,Wi,1000);
    if extend_velocities==1
        G=vel>-400 & vel<400;
    else
        G=vel>-200 & vel<200;
    end
    U=vel(G);
    LSDy=LSD(:,G);
    coarse_ccf.velocity=U;
    coarse_ccf.intensity=LSDy;
    coarse_ccf.jd=jd;
    eval(['save(''ccf_' objectname '.mat'',''coarse_ccf'')'])
    if supress_figures==0
        figure
        plot(U,LSDy)
        eval(['title(''Line Profiles ' objectname(4:end) ''')'])
    end
end


%OPTION:Radial Velocity & vsini
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do a coarse CCF to get a systemic velocity
if cross_correlation==1 && rv_measurement==1
    fprintf('Measuring radial velocity and vsini\n')
    [systemic_velocity_direct,systemic_velocity_gauss]=radial_velocity_calculation_systemic_only(coarse_ccf);
    finaldata.systemic_velocity=systemic_velocity_direct;
    finaldata.systemic_velocity_gauss=systemic_velocity_gauss;
end

%OPTION:Cross-correlation Part 2 - least-squares fit delta function (needs systemic velocity measured above)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cross_correlation==1
    %FINALISING CCFS FROM PART 1
    %getting weights from a noise estimate and producing a summed spectrum
    intfilt=medfilt1(fullint(:,1:10000)',11)';
    resid=fullint(:,1:10000)-intfilt;
    weight=1./std(resid,[],2);
    
    %LEAST SQUARED DELTA FUNCTION CCF
    %remove very weak lines
    SS2=Wi<weak_limit2;
    W(SS2)=[];
    Wi(SS2)=[];
    delfn=cat(2,W(:),Wi(:));
    
    %remove the systemic velocity from the spectrum before fitting
    if exist('systemic_velocity_direct','var') == 0 || isnan(systemic_velocity_direct) || rv_measurement==0
        systemic_velocity_direct=0;
    end
    [newwave]=shift_velocity(fullwave,systemic_velocity_direct);
    %create least-squares matched delta function and CCF
    if size(fullint,1)==1
        fullint_mean=fullint;
    else
        fullint_mean=mean(fullint);
    end
    newdelfn = LS_match_delfn(newwave,fullint_mean,delfn,65,27.3,0.35);% default 65,27.3,0.35
    [v2,vel2,LSD2,CCF2,dfn_auto2,wcw2] = do_ccf_from_spec_wavelengths_and_widths(newwave,fullint,newdelfn(:,1),newdelfn(:,2),1000);
    
    %plotting
    if extend_velocities==1
        G2=vel>-400 & vel<400;
    else
        G2=vel>-200 & vel<200;
    end
    U2=vel2(G);
    LSDy2=LSD2(:,G);
    
    %plot the LSDs
    if supress_figures==0
        figure
        plot(U2,LSDy2)
        eval(['title(''Least squares fit CCF ' objectname(4:end) ''')'])
    end
    least_squares_fit_ccf.velocity=U2;
    least_squares_fit_ccf.intensity=LSDy2;
    least_squares_fit_ccf.jd=jd;
    eval(['save(''least_squares_fit_ccf_' objectname '.mat'',''least_squares_fit_ccf'')'])
    LSDyy2=LSDy2+(1-max(mean(LSDy2)));
    
    %RVS from LEAST SQUARED CCF
    if rv_measurement==1
        [radial_velocity2,radial_velocity_error2,systemic_velocity_direct2,systemic_velocity_gauss2,vsini2,vsini_error2]=radial_velocity_calculation(least_squares_fit_ccf);
        
        supplementary_data.radial_velocity_lsf=radial_velocity2;
        supplementary_data.radial_velocity_error_lsf=radial_velocity_error2;
        supplementary_data.systemic_velocity_direct_lsf=systemic_velocity_direct2;
        supplementary_data.systemic_velocity_gauss_lsf=systemic_velocity_gauss2;
        supplementary_data.vsini=vsini2;
        supplementary_data.vsini_error=vsini_error2(1);
        
        %calculate single vsini
        supplementary_data.single_vsini=mean(vsini2);
        supplementary_data.single_vsini_error=vsini_error2(1)/sqrt(numel(vsini2));
        supplementary_data.single_vsini_std=std(vsini2);
        eval(['save(''supplementary_results_' objectname '.mat'',''supplementary_data'')'])
    end
    %RVS from standard CCF
    if rv_measurement==1
        [radial_velocity3,radial_velocity_error3,systemic_velocity_direct3,systemic_velocity_gauss3,vsini3,vsini_error3]=radial_velocity_calculation(coarse_ccf);
        
        finaldata.radial_velocity=radial_velocity3;
        finaldata.radial_velocity_error=radial_velocity_error3;
        finaldata.vsini=vsini3;
        finaldata.vsini_error=vsini_error3(1);
        
        %calculate single vsini
        finaldata.single_vsini=mean(vsini3);
        finaldata.single_vsini_error=vsini_error3(1)/sqrt(numel(vsini3));
        finaldata.single_vsini_std=std(vsini3);
    end
    
end

% try
%     finaldata = rmfield(finaldata,["vsini_direct", "vsini_gauss", "vsini_gauss_error", "single_vsini_lsf", "single_vsini_error_lsf", "single_vsini_std_lsf", "vsini_gauss_lsf","vsini_gauss_error_lsf","vsini_error_lsf","vsini_direct_lsf" ]);
% catch
% end

%FINALISING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order_merge==1
    if supress_figures==0
        figure
        plot(fullwave,fullint)
        eval(['title(''Full Spectra ' objectname(4:end) ''')'])
    end
else
    if supress_figures==0
        figure
        plot(finaldata.wave87,finaldata.int87)
        eval(['title(''Sample Spectrum, order 87, ' objectname(4:end) ''')'])
    end
    
end

if rv_measurement==1 && all(isnan(radial_velocity3)~=1)
    if supress_figures==0
        figure
        errorbar(jd,radial_velocity3,radial_velocity_error3*ones(numel(jd),1),'xb')
        if isnan(radial_velocity2)~=1
            hold on
            errorbar(jd,radial_velocity3,radial_velocity_error3*ones(numel(jd),1),'xr')
        end
        eval(['title(''Radial Velocity ' objectname(4:end) ''')'])
    end
else
end

if rv_measurement==1
    if all(isnan(vsini2)~=1)
        if supress_figures==0
            figure
            errorbar(jd,vsini3,vsini_error3,'x')
            if isnan(vsini2)~=1
                errorbar(jd,vsini3,vsini_error3,'x')
            end
            eval(['title(''vsini ' objectname(4:end) ''')'])
        end
    else
    end
end

if FAMIAS_output==1 && cross_correlation==1
    fprintf('Writing files for FAMIAS\n')
    warning off
    mkdir('CCF_FAMIAS_files')
    cd('CCF_FAMIAS_files')
    eval(['write_files_for_FAMIAS(jd,U,LSDy,weight,''times_' objectname(4:end) '.txt'',''' objectname(4:end) '_'')'])
    cd ../
    mkdir('LSF_CCF_FAMIAS_files')
    cd('LSF_CCF_FAMIAS_files')
    eval(['write_files_for_FAMIAS(jd,U2,LSDy2,weight,''times_' objectname(4:end) '.txt'',''' objectname(4:end) '_'')'])
    cd ../
    fprintf('FAMIAS-ready files written in reduced data directory\n')
    warning on
end

myflag=1;

finaldata.signal_to_noise=signal_to_noise;

if exist('finaldata','var')==0
    finaldata=[];
end

fprintf('DATA PROCESSED\n')

