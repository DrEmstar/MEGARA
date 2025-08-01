function psynthint=get_synth_piece(synthwave,synthint,wave)
minwave=min(wave);maxwave=max(wave);step=round(mean(diff(wave))*1000)/1000;
inds=findnearest(synthwave,minwave-5*step);
inde=findnearest(synthwave,maxwave+5*step);
tsynthw=synthwave(inds:inde);
tsynthi=synthint(inds:inde);
psynthint=spline(tsynthw,tsynthi,wave);
