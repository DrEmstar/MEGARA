function bcorr=get_barycentric_correction(starname,expdate,midtime,filename)
z=pwd;
pat=which('get_barycentric_correction.m');
eval(['cd ' pat(1:end-28)]); 
setenv('HRSP','./')
bcorr=evalc(['!./barycorr_hrsp ' starname ' ' expdate ' ' midtime]);
if size(bcorr,2) == 9
else
    header = get_header(filename);
    header2=get_required_headers_from_header(header);
    fid=fopen('stars.dat','w');
    try
        fprintf(fid,'1 ICRS 2000.0 %7s %7s 6.9514 -15.895 17.87 50.799 \n',header2.RA_2000,header2.DEC_2000);
    catch
        fprintf(fid,'1 ICRS 2000.0 %7s %7s 6.9514 -15.895 17.87 50.799 \n',header2.RA,header2.DEC);
    end
    fclose(fid);
    bcorr=evalc(['!./barycorr_hrsp USR_1 ' expdate ' ' midtime]);
    if size(bcorr,2) == 9
    else
    warning('Barycentric correction not calculated correctly!')
    end
end
cd(z)