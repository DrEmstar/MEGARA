function [final_delfns, residual, thefit] = leastsquaresfit_spec_w_broad_delfns(wavelength,spec,delfns,vsini,intrinFWHM,epsx)

%everythings in columns
wavelength=wavelength(:);
spec=spec(:);

%delfns must be a 2-column list of wavelength and equivalent width in mA

delfnwavs=delfns(:,1);

centralwave=mean(wavelength);
step = round(mean(diff(wavelength))*1000)/1000;

%we have: vsini, intrinFWHM, epsx
result = make_spectral_line(vsini,intrinFWHM,epsx,centralwave,step);
%make line have equivalent width of 1 and invert
r=1-result(:,2);
r=r/(sum(r*step));

save('temp_delfn_params.mat','r','delfnwavs')

x0=delfns(:,2);
lb = zeros(size(x0));
ub = ones(size(x0))*500;

options=optimset('Display','off');
%do the least-squares fitting
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@delfn_fitting,x0,wavelength,spec,lb,ub,options);
final_delfns=cat(2,delfnwavs(:),x(:));
thefit = delfn_fitting(x,wavelength);

function result = make_spectral_line(vsini,intrinFWHM,epsx,centralwave,step)
%input parameters are:
%vsini        vsini of synth line in km/s
%intrinFWHM   intrinsic FWHM of synth line in km/s
%epsx         limb-darkening parameter
%centralwave  central wavelength used to obtain profile being measured
%outputwave   wavelength vector to output (project) the synth line on to

c=299792.458; %speed of light (km/s)
% Define rotational broadening curve over fine wavelength grid initially
deltalamda=-10:0.001:10;
%prepare parameters for broadening curve
lamda=centralwave;
lam=deltalamda/lamda;lamgrp=lam/(vsini/c);
%rotational broadening fn from Gray 1992
G=real((2*(1-epsx)*((1-lamgrp.^2).^0.5)+0.5*pi*epsx*(1-lamgrp.^2))/(pi*lamda*vsini/c*(1-epsx/3)));
indices=find(G > 0);  %indices to the positive part of G
if numel(indices)/2 == round(numel(indices)/2)  %if even number of elements left in G
    %must remove the one closest to zero as we want an odd number of elements
    ttt=find([G(indices(1)) G(indices(end))] ==min([G(indices(1)) G(indices(end))]));
    if ttt==1
        indices(1)=[];
    else
        indices(end)=[];
    end
end
broadcurve=G(min(indices):max(indices));%select the right parts of G
% make gaussian curve to be broadened in same wavelength space as broadening profile
a1=1; b1=0;
%convert FWHM (velocity to wavelength)
FWHM=intrinFWHM; FWHM=FWHM*lamda/c;
%make the gaussian
gaussiancurve=a1*exp(-(2*sqrt(log(2))*(deltalamda-b1)/FWHM).^2);
% convolve the gaussian line and the broadening curve to produce the
% broadened curve
res=real(conv(broadcurve,gaussiancurve));
% scale to height 1 and invert
res=res/max(res); res=1-res;
% make x axis for result
%ceil(numel(res)/2); is the central element
newax=1:numel(res); %create
newax=newax-ceil(numel(res)/2); %centre
newax=newax*0.001; %scale
%trim to just at line limits
d=find(res<0.999);
res=res(min(d)-100:max(d)+100);     %100 pixels is 0.1A
newax=newax(min(d)-100:max(d)+100);
%make axis based on step input
st=ceil(newax(1)/step)*step;
ed=floor(newax(end)/step)*step;
finax=st:step:ed;
%project on to requested wavelength axis
newres=spline(newax,res,finax);
%combine xax and int vector
result=finax(:);
result(:,2)=newres(:);
warning off
delete temp_delfn_params.mat
warning on
