function bcorr=get_barycentric_correction(starname,expdate,midtime,filename)
z=pwd;
pat=which('get_barycentric_correction.m');
eval(['cd ' pat(1:end-28)]); 
bcorr=evalc(['!./barycorr_HRSP ' starname ' ' expdate ' ' midtime]);
if size(bcorr,2) == 9
else
    disp('WARNING: Barycentric correction not calculated correctly!')
end
cd(z)