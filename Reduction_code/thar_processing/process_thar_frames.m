function [jd, m, fwave, numbad_thar_flag]=process_thar_frames(tharname, info, allorders, Bias, flatdata, blue_data_chop_value)

%ThAr image processing (1 image at a time only)
backfit.back_polyfit = 0; %ThAr images don't have background measurements
backfit.xscale=1;backfit.xcen=0;
backfit.yscale=1;backfit.ycen=0;
backfit.zscale=1;backfit.zcen=0;

sname=tharname(1:end-4);
eval(['jd=info.' sname '.MJD_OBS;'])
data = raw_hercules_fitsread(tharname);

%assuming no Bias image
if ~isempty(Bias)
    data=data-Bias;
end
data = blue_data_chop(data,blue_data_chop_value);
thardata= extract_all_orders_no_background(allorders,data);
for s=1:flatdata.numords
    eval(['thardata.order_' num2str(s) '.summed_data=sum(thardata.order_' num2str(s) '.data(:,flatdata.order_' num2str(s) '.extraction_inds),2);'])
end

%Handling rollover of 4 digit filenames from J5XXXX to J6XXXX

if jd<60000

%Dynamic thardat selection tested for 2007-2018 data

thardate=str2num(sname(2:5));

if thardate <= 4463
    load thardat_g.mat
elseif (thardate > 4463) && (thardate <= 4781)
    load thardat_f.mat
elseif (thardate > 4781) && (thardate <= 5454)
    load thardat_d.mat
elseif (thardate > 5454) && (thardate <= 5702)
    load thardat_e.mat
elseif (thardate > 5702) && (thardate <= 6337)
    load thardat_c.mat
elseif (thardate > 6337) && (thardate < 6791)
    load thardat_b.mat
elseif (thardate > 6791) && (thardate < 7083)
    load thardat_a.mat
elseif (thardate >= 7083) && (thardate < 7178)
    load thardat_i.mat
elseif (thardate >= 7178) && (thardate < 7821)
    load thardat_a.mat
elseif (thardate >= 7821) && (thardate < 7968)
    load thardat_h.mat
elseif (thardate >= 7968) && (thardate < 8120)
    load thardat_a.mat
elseif thardate >= 8120 && (thardate < 8144)
    load thardat_j.mat
elseif thardate >= 8144 && (thardate < 8146)
    load thardat_k.mat
elseif thardate >= 8146 && (thardate < 8303)
    load thardat_j.mat
elseif thardate >= 8303 && (thardate < 9461)
    load thardat_l.mat
elseif thardate >= 9461 && (thardate < 9650)
    load thardat_mm.mat
   elseif thardate >= 9650
       load thardat_n.mat
end

else 
    load thardat_o.mat
end



try
    isempty(air);
catch
    try
        fid=fopen('thar.dat','r');
        C=textscan(fid,'%.0f %.4f %.4f %.4f %s %.3f 2017%.3f %.3f %.3f %.3f %.3f');
        order=C{1};air=C{2};vac=C{3};wavnum=C{4};species=C{5};ul=C{6};uc=C{7};ur=C{8};vl=C{9};vc=C{10};vr=C{11};
        status=fclose(fid);clear status
        % standard corrections to the text file
        uc=uc/0.015+36;vc=vc/0.015-200;
        ul=ul/0.015+36;vl=vl/0.015-200;
        ur=ur/0.015+36;vr=vr/0.015-200;
    catch
        disp('there is no thardat.mat or thar.dat file accessible')
        disp('without these calibrations ThAr processing cannot continue')
        return
    end
end

%apply simple program to change to absolute order numbers
[thardata_2 m]=redo_order_numbering_simple(thardata,allorders);

%finding and fitting ThAr lines in extracted spectrum
%make plot_q=1 if you want plots of every th line as it is processed to
%display on screen pausing for the enter button to continue
plot_q=0;
[thfitinfo,numbad_thar_flag] = find_th_lines(thardata_2,order,ul,uc,ur,plot_q);

%wavelength calibration
wavfit=wavelength_calibration(thfitinfo,order,air);
%make wavelength axes for each ThAr for data:
fwave=make_wave_vector(thardata_2,wavfit,m);

end