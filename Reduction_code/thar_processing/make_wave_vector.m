function wave=make_wave_vector(thardata,wavfit,m)

for s=1:numel(m)
    eval(['xpos=thardata.order_' num2str(m(s)) '.xax;']); xpos=xpos(:);
    yy=(xpos-wavfit.ycen)/wavfit.yscale;
    xx=(m(s)-wavfit.xcen)/wavfit.xscale;
    twav=polynomial_val_2d(wavfit.p,repmat(xx,size(yy)),yy);
    twav=twav*wavfit.zscale+wavfit.zcen;
    eval(['wave.order_' num2str(m(s)) '=twav/m(s);'])
end
