function [datacf] = manual_continuum_fitting(teff,logg,vsini,sigm,wstart,wstop,shift,firstorder,lastorder,data2,objectname)
%This program allows manual identification of the continuum level for each
%order. Try to keep the function smooth and fit the top sides of the line
%profiles.

[synthwave, synthint, wavelengths, ele, elem, widths] = make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop);
synthwave=shift_velocity(synthwave,shift);
for ord=[firstorder:lastorder]
    eval(['int=data2.int' num2str(ord) ';'])
    eval(['wave=data2.wave' num2str(ord) ';'])
    temp_synthint=spline(synthwave,synthint,wave);
    if length(median(int))==1
        med_int=int;
    else
        med_int=median(int);
    end
    [result,thecfit] = continuum_fitting_inline(wave,med_int,temp_synthint);
    eval(['datacf.int' num2str(ord) '=thecfit;'])
    eval(['datacf.w' num2str(ord) '=wave;'])
end

customdatacf=eval(['''datacf_' objectname '''']);
save(customdatacf,'datacf')
