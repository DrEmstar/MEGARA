function [synthwave, synthint, wavelengths, ele, elem, widths] = make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop)
z=pwd;
%change to appropriate synspec path
pat=which('synspec_pointer.m');
eval(['cd ' pat(1:end-17)]); 
%round Teff to nearest 500K & logg to nearest 0.5
teff_500=round(teff/500)*500;
logg_05=round(logg/0.5)*0.5;
write_new_file_for_synspec(teff_500,logg_05);
write_mod55_file_for_synspec(wstart,wstop);
evalc('!./rmod');
evalc('!./KSynspec output mod gf3000');
[wavelengths, ele, elem, widths] = read_in_output12;
write_r_file_for_rotin3(vsini,sigm,wstart,wstop);
evalc('!./rotin3 < r.dat');
%check everything is correct and writing new files
syninfo=dir('output.11');
file_time=datetime(syninfo.date);
now_time=datetime('now');
test_time=now_time-minutes(5);
if file_time < test_time
    warning('SynSpec reading previous spectrum. Check star parameters')
else
end
[synthwave, synthint] = read_in_output11;
cd(z)