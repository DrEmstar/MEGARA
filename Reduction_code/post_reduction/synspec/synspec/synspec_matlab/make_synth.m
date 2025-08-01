function [synthwave, synthint, wavelengths, ele, elem, widths] = make_synth(teff,logg,vsini,sigm,wstart,wstop)
z=pwd;
%change to appropriate synspec path
pat=which('synspec_pointer.m');
eval(['cd ' pat(1:end-17)]); 
write_new_file_for_synspec(teff,logg);
write_mod55_file_for_synspec(wstart,wstop);
evalc('!./rmod');
evalc('!./KSynspec output mod gf3000');
[wavelengths, ele, elem, widths] = read_in_output12;
write_r_file_for_rotin3(vsini,sigm,wstart,wstop);
evalc('!./rotin3 < r.dat');
[synthwave, synthint] = read_in_output11;
cd(z)
