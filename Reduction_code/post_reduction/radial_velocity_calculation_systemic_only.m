function [systemic_velocity_direct,systemic_velocity_gauss]=radial_velocity_calculation_systemic_only(coarse_ccf)

%determine extent of profile
summed_line=sum(coarse_ccf.intensity,1);
[pks, locs, w]=findpeaks(max(summed_line)-summed_line,coarse_ccf.velocity,'Annotate','extents','MinPeakProminence',max(max(summed_line)-summed_line)/2,'WidthReference','halfheight');
[negpks, neglocs]= findpeaks(summed_line,coarse_ccf.velocity);
try
minidx = find(neglocs < locs-w/1.2, 1, 'last');
catch
    systemic_velocity_direct=NaN;
    systemic_velocity_gauss=NaN;
    return
end
if isempty(minidx)
    minloc=neglocs(1);
else
    minloc = neglocs(minidx);
end
maxidx = find(neglocs > locs+w/1.2, 1, 'first');
if isempty(maxidx)
    maxloc=neglocs(end);
else 
    maxloc = neglocs(maxidx);
end

for n=1:length(coarse_ccf.jd)
    botpix=find_nearest(coarse_ccf.velocity,minloc);
    toppix=find_nearest(coarse_ccf.velocity,maxloc);
    bkgrd=mean([coarse_ccf.intensity(n,1:(botpix-10)) coarse_ccf.intensity(n,(toppix+10):end)]);
    try
    [vsini(n),fittedint]=vsini_directfit(coarse_ccf.velocity,coarse_ccf.intensity(n,:)./bkgrd,locs);
    [v_gpks(n),v_glocs(n),v_gwidths(n),v_gproms(n),bPk,iLB,iRB]=findpeaks_vsini(1-fittedint,coarse_ccf.velocity,'Annotate','extents','MinPeakProminence',max(1-fittedint)/10,'WidthReference','halfprom');
    catch
        v_glocs(n)=NaN;
        continue
    end
end

systemic_velocity_direct=locs;
systemic_velocity_gauss=nanmedian(v_glocs);