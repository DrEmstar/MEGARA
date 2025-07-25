function V = vsini_fittingfn(x,xdata);
%works out broadened line profile in lamda space and converts to velocity space
load temp_epsx.mat
load temp_lambda0.mat
warning off
% x(1)=FWHM gaussian
% x(2)=V sin(i)
c=299792.458; %speed of light (km/s)
deltalamda=-10:0.001:10; %make a fine wavelength grid on which to construct
deltalamda1=-10:0.001:10; %our broadening curve
indices=find(real((2*(1-epsx)*((1-((deltalamda*c/lamda0)/(x(2))) .^2).^0.5)+0.5*pi*epsx*(1-((deltalamda*c/lamda0)/(x(2))) .^2)) / (pi*lamda0*x(2)/c*(1-epsx/3)) ) > 0);
start=indices(1);
finish=indices(end);
deltalamda=deltalamda(start:finish);
K=real(conv(real((2*(1-epsx)*((1-((deltalamda*c/lamda0)/(x(2))).^2).^0.5)+0.5*pi*epsx*(1-((deltalamda*c/lamda0)/(x(2))).^2))/(pi*lamda0*x(2)/c*(1-epsx/3))),exp(-(2*sqrt(log(2))*deltalamda1/x(1)).^2)));
centre1=find(deltalamda==0);
centre2=find(deltalamda1==0);
centre=centre1+centre2-1;
%convert width parameter to velocity units
%make an x-axis for the convolution
newxax=-(centre-1)*0.001:0.001:0;
t=numel(K)-numel(newxax);
newxax=-(centre-1)*0.001:0.001:t*0.001;
%convert the x-axis to velocity space
newxax=newxax/lamda0*c;
%project it on to the required axis
V=spline(newxax,K,xdata);
V=V/max(V); %normalise
% hold on
% plot(newax,V)
warning on