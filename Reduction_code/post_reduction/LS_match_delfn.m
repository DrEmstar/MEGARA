function [newdelfn,thefits,waves] = LS_match_delfn(wavelength,spec,delfn,vsini,intrinFWHM,epsx)

start = delfn(1,1);
cnt=0;
crossover = 7; %size of crossover wavelength in each piece fitted (7 is good enough)
stepsize = 50; %size of wavelength region to be fitted at one time (50 is pretty good with reasonable fitting times)
while start < delfn(end,1)
    finish = start+stepsize;
    delind = find(delfn >= start & delfn < finish);
    if isempty(delind)
        %no lines in this region
        start=finish;
        continue
    end
    cnt=cnt+1;
    tdelfn = delfn(delind,:);
    waveind = find(wavelength > start-1 & wavelength < finish+1);
    wave = wavelength(waveind);sp=spec(waveind);
    [final_delfns, residual, thefit] = leastsquaresfit_spec_w_broad_delfns(wave,sp,tdelfn,vsini,intrinFWHM,epsx);
    eval(['thefits.f' num2str(cnt) '=thefit;'])
    eval(['waves.w' num2str(cnt) '=wave;'])
    if cnt == 1
        %chop end (not start) off since there is overlap
        wantedinds=find(final_delfns(:,1) < finish-(crossover/2));
        final_delfns=final_delfns(wantedinds,:);
        newdelfn=final_delfns;
    else
        %chop both ends off to get rid of poor end fitting (this is why there is
        %overlap between fits)
        wantedinds=find(final_delfns(:,1) >= start+(crossover/2) & final_delfns(:,1) < finish-(crossover/2));
        final_delfns=final_delfns(wantedinds,:);
        newdelfn=cat(1,newdelfn,final_delfns);
    end
    start = finish-crossover;
end

