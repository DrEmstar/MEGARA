function stellar_frame_processing(stellarimgname,info,allorders,flatdata,fwave1,fwave2,jd1,jd2,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,jdexp)
%STELLAR image processing
sname=stellarimgname(1:end-4);
try
    eval(['jdexp=info.' sname '.MJD_OBS;'])
catch
    disp('couldn''t get Julian date, check filename and')
    disp('info file and try again')
    return
end
data = raw_hercules_fitsread(stellarimgname);
data = blue_data_chop(data,blue_data_chop_value);

warning off
backfit= fit_background(data,allorders);
warning on

stardata = extract_all_orders_no_background(allorders,data);
backdata = extract_all_orders_no_background(allorders,backfit);

%cosmic ray filtering
stardata.numords=allorders.numords;
stardata = remove_cosmics(stardata,cosmic_plotting);

for s=1:flatdata.numords
    eval(['backdata.order_' num2str(s) '.summed_data=sum(backdata.order_' num2str(s) '.data(:,flatdata.order_' num2str(s) '.extraction_inds),2);'])
    eval(['backdata.order_' num2str(s) '.summed_smoothed_data=mean_smoothing(backdata.order_' num2str(s) '.summed_data,backdata_medfilt);'])
    
    eval(['stardata.order_' num2str(s) '.summed_data=sum(stardata.order_' num2str(s) '.data(:,flatdata.order_' num2str(s) '.extraction_inds),2);'])
    eval(['stardata.order_' num2str(s) '.summed_data_b=stardata.order_' num2str(s) '.summed_data-backdata.order_' num2str(s) '.summed_smoothed_data;'])
end

%apply simple program to change to absolute order numbers
[stardata m]=redo_order_numbering_simple(stardata,allorders);

%make wavelength axes for stellar obs by linear interpolation of ThAr axes
wavelength = linear_interpolate_wavelengths1D(fwave1,fwave2,jd1,jd2,jdexp,m);
for s=1:numel(m)
    eval(['fl=flatdata.order_' num2str(s) '.sumdata;'])
    flb=fl;
    fl=medfilt1(fl,flatdata_medfilt); %median filter the flat-field data
    fl=fl/median(fl);F=fl==0;fl(F)=1e-6; %normalise by max value and set zeros to low value

    eval([sname '.order_' num2str(m(s)) '(:,1)=wavelength.order_' num2str(m(s)) ';'])
    eval(['xxx=stardata.order_' num2str(m(s)) '.summed_data_b(:);'])
    xxx=xxx/median(xxx);
    xxxf=xxx(:)./fl(:);
    eval([sname '.order_' num2str(m(s)) '(:,2)=xxxf;'])
    
    ind1=find_nearest(fl(1:round(numel(fl)/2)),0.1);
    ind2=find_nearest(fl(round(numel(fl)/2):end),0.1);
    ind2=ind2+round(numel(fl)/2)-1;
    
    eval([sname '.order_' num2str(m(s)) '=' sname '.order_' num2str(m(s)) '(ind1:ind2,:);'])
    
    eval([sname '_prc.order_cutpixels.order_' num2str(m(s)) '.ind1=ind1;'])
    eval([sname '_prc.order_cutpixels.order_' num2str(m(s)) '.ind2=ind2;'])
    
    eval([sname '_prc.order_bkgd_' num2str(m(s)) '(:,1)=backdata.order_' num2str(s) '.summed_data;'])
    eval([sname '_prc.order_bkgd_smoothed_' num2str(m(s)) '(:,1)=backdata.order_' num2str(s) '.summed_smoothed_data;'])
    
    eval([sname '_prc.order_flat_' num2str(m(s)) '(:,1)=flb;'])
    eval([sname '_prc.order_flat_smoothed_' num2str(m(s)) '(:,1)=fl;'])
    
    eval([sname '_prc.order_nff_' num2str(m(s)) '(:,1)=stardata.order_' num2str(m(s)) '.summed_data;'])
    eval([sname '_prc.notes=''Reduced spectrum is: ($file_prc.order_nff_ - $file_prc.order_bkgd_smoothed_) / ($file_prc.order_flat_smoothed_). background smoothing is at medfilt =' num2str(backdata_medfilt) ', flat smoothing is at ' num2str(flatdata_medfilt) ''';'])
end

eval([sname '.jd=jdexp;'])
eval(['header=get_header(''' sname '.fit'');'])
headers = get_required_headers_from_header(header);
starname=headers.OBJECT;
expdate=headers.DATE_OBS;  % OLD CCD SOFTWARE
if expdate(12)=='T'
    expdate=headers.DATE;% NEW CCD SOFTWARE
end
starttime=headers.REC_STRT;
if str2num(starttime(end-4:end-1))>=60
    starttime(end-4:end-1)='59.9';
end
  
%calculating signal-to-noise
stellar_order=stardata.order_100;
exp_time=headers.EXPTIME;
signal_to_noise=compute_signal_to_noise_spectrum(stellar_order,exp_time);
eval([sname '.signal_to_noise=signal_to_noise;'])

midexp=headers.HERCFWMT;
if isempty(midexp) || isnan(midexp)
    warning on
    warning('Flux-weighted mid-time unavailable for %s, using exposure mid-time instead',sname)
    midexp=exp_time/2;
end
%Adding the mid exposure time
t1=datetime(expdate(2:end-1),'InputFormat','yyyy-MM-dd');
t2=datetime(starttime(2:end-1),'InputFormat','HH:mm:ss.S');
tt=t2;tt.Year=t1.Year;tt.Month=t1.Month;tt.Day=t1.Day;
midtime=tt+seconds(midexp);
%putting in Terrestial Time instead of UTC for hrsp_barycorr
midtime_tt=midtime+seconds(69.184);
midtime_utc=['"' datestr(midtime,'HH:MM:SS.FFF') '"'];
midtime_tt=['"' datestr(midtime_tt,'HH:MM:SS.FFF') '"'];

%saves some headers for adding to the '.mats'
fibre = headers.HERCFIB;     %Hercules Fibre 1/2/3
counts = headers.HERCFTC;     %Hercules flux meter total counts

try
    RA_j2000=headers.RA; %RA in j2000
catch
    RA_j2000=headers.RA_2000;
end
if isempty(RA_j2000)==1
    try 
        RA_j2000=headers.RA_2000;
    catch
        warning('RA information missing')
    end    
end

try
    DEC_j2000=headers.DEC; %DEC in j2000
catch
    DEC_j2000=headers.DEC_2000;
end

if isempty(DEC_j2000)==1
    try 
        DEC_j2000=headers.DEC_2000;
    catch
        warning('DEC information missing')
    end
end


eval([sname '.starname=starname;'])
eval([sname '.expdate=expdate;'])
eval([sname '.midtime_utc=midtime_utc;'])
eval([sname '.midtime_tt=midtime_tt;'])
eval([sname '.sname=sname;'])
eval([sname '.blue_data_chop=blue_data_chop_value;'])
eval([sname '.counts=counts;'])
eval([sname '.fibre=fibre;'])
eval([sname '.exptime=exp_time;'])
eval([sname '.blue_data_chop=blue_data_chop_value;'])
eval([sname '.RA_j2000=RA_j2000;'])
eval([sname '.DEC_j2000=DEC_j2000;'])

fprintf('saving %s.mat and %s.mat\n\n',sname, [sname '_prc'])
save(sname,sname)
save([sname '_prc'],[sname '_prc'])
% Get list of variable names in current function workspace
vars = who;
% Clear them all 
clearvars(vars{:});