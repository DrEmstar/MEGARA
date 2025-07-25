function write_mod55_file_for_synspec(wstart,wstop)
z=pwd;
%change to appropriate synspec path
pat=which('synspec_pointer.m');
eval(['cd ' pat(1:end-17)]); 
fid = fopen('mod.55','w+');

fprintf(fid,'       0      55       0                    ! imode,idstd,iprint\n');
fprintf(fid,'       0       0       0       0            ! imodel,ichang,interp,ichemc\n');
fprintf(fid,'       0                                    ! Lyman quasi-mol\n');
fprintf(fid,'       0       0       0       0       0\n');
fprintf(fid,'     -23      24      26                    ! line broadening (H,HeI,HeII)\n');
fprintf(fid,'    %.1f   %.1f    10       0  1.d-3    0.1\n',wstart,wstop);
fprintf(fid,'      0\n');
fprintf(fid,'      2.0\n');
fclose(fid);
cd(z)
