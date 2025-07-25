function [centre, FWHM] = find_max_and_FWHM(xdata, ydata)
ydata=ydata-min(ydata);
ydata=ydata/max(ydata);
indcen=xdata(round(numel(xdata)/2));

[pks,locs,widths,proms]=findpeaks(ydata,xdata,'MinPeakHeight',0.4,'WidthReference','halfheight');
if numel(pks)>1
    centre=locs(find_nearest(locs,indcen));
    FWHM=widths(find_nearest(locs,indcen));
elseif numel(pks)==1
    centre=locs;
    FWHM=widths;
else
    centre=0;
    FWHM=0;
end
