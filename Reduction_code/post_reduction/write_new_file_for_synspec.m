function write_new_file_for_synspec(teff,logg)
%change to appropriate synspec path
pat=which('synspec_pointer.m');
eval(['cd ' pat(1:end-17)]); 

fid = fopen('new','w+');
if teff < 10000
    fprintf(fid,'       %.0f.          %.2f',teff,logg);
else
    fprintf(fid,'      %.0f.          %.2f',teff,logg);
end
fclose(fid);
