function sspec = delfn_fitting(x,xdata)

%x is the list of delfn depths
%xdata is the 
%need to have also the delfn wavelengths
%global delfnwavs r
load temp_delfn_params

rnddelfnwave=round(delfnwavs*10)/10;
[b, m, n] = unique(rnddelfnwave);
delfn_array=zeros(size(xdata));
for s=1:numel(b)
    delfn_array(xdata==b(s))=sum(x(n==s));
end

%convolve with delfn
S=1-conv(delfn_array/1000,r); %delfn/1000 since they are in milli-Angstroms and we're in Angstroms
%trim to input size
sspec=S(floor(numel(r)/2)+1:end-floor(numel(r)/2));
