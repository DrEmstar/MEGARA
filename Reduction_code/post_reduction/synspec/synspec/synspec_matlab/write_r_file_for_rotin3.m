function write_r_file_for_rotin3(vsini,sigm,wstart,wstop)
z=pwd;
%change to appropriate synspec path
pat=which('synspec_pointer.m');
eval(['cd ' pat(1:end-17)]); 
fid = fopen('r.dat','w+');

fprintf(fid,' ''output.7''   ''output.17''    ''output.11'' \n');
fprintf(fid,'     %.0f.    0.05       0.1\n',vsini);
fprintf(fid,'      %.2f             0\n',sigm);
fprintf(fid,'    %.0f.   %.0f.      1\n',wstart,wstop);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fclose(fid);
cd(z)
