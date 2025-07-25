function [v,vel,LSD,CCF,dfn_auto,wcw] = do_ccf_from_spec_wavelengths_and_widths(wave,intmtx,wa,wi,numbins)

c = 299792.458; %speed of light (km/s)
wcw=sum(wi/sum(wi).*wa);

%go to log space
logwave=log(wave);
minlogaxis=floor(min(logwave)*100000)/100000;
maxlogaxis=ceil(max(logwave)*100000)/100000;
logaxisstep=(maxlogaxis-minlogaxis)/length(wave);
logaxis=minlogaxis:logaxisstep:maxlogaxis;
velstep=(exp(log(wcw)+logaxisstep)-wcw)/wcw*c;
lwa=log(wa);
ldelfn = cat(2,lwa(:),wi(:));
ldelfn_array = logdelfn_2_logdelfn_array(ldelfn,logaxis);
[dfn_auto,xax] = xcorr(ldelfn_array,numbins,'unbiased');
v=xax*velstep;
vel=ceil(min(v)):velstep:floor(max(v));

for s=1:numel(intmtx(:,1))
    % make spec linear in log
    logspecint=spline(logwave,intmtx(s,:),logaxis);
    CCF(s,:) = xcorr(1-logspecint,ldelfn_array,numbins,'unbiased');
    tresult = testing_deconvolution(CCF(s,:),dfn_auto);
    tresult = (1-tresult)/(max(1-tresult));
    LSD(s,:) = spline(v,tresult,vel);
end

