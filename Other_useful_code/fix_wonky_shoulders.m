%Put star information here
clear all
objectname=['hd_175362'];


load ("final_data_" + objectname + ".mat")
load ("ccf_" + objectname + ".mat")
fullwave=finaldata.fullwave;
fullint=finaldata.fullint;
%jd=finaldata.jd;


figure
plot(coarse_ccf.velocity,coarse_ccf.intensity)
[x,y] = getpts;
p = polyfit(x,y,1);
f1 = polyval(p,coarse_ccf.velocity);
norm_int=coarse_ccf.intensity./f1;
hold on
plot(coarse_ccf.velocity,norm_int)
coarse_ccf.intensity=norm_int;


%getting weights from a noise estimate and producing a summed spectrum

%Find Nans and set to 1
G=isnan(fullint);
fullint(G)=1;

intfilt=medfilt1(fullint(:,1:10000)',11)';
resid=fullint(:,1:10000)-intfilt;
weight=1./std(resid,0,2);

U=coarse_ccf.velocity;
LSDy=coarse_ccf.intensity;
jd=coarse_ccf.jd;


cd("Reduced_Data/hd_" + objectname(4:end) + "/")
mkdir("Wonky_fix_CCF_FAMIAS_files_" + objectname(4:end))
cd("Wonky_fix_CCF_FAMIAS_files_" + objectname(4:end))
eval(['save(''ccf_' objectname '_wonky_fix.mat'',''coarse_ccf'')'])
eval(['write_files_for_FAMIAS(jd,U,LSDy,weight,''times_' objectname(4:end) '.txt'',''' objectname(4:end) '_'')'])
