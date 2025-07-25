function [vsini,fittedint,int]=vsini_directfit(vel,int,locs);
%format long g
%global lamda0 epsx
% Fit the line profile by optimizing an unconstrained FWHM Gaussian convolved with the rotational broadening curve of Gray 1992 for the specified wavelength
% (usually the weighted centre wavelength of the ccf-function). Note that the input spectrum is in velocity space.
c=299792.458; % Speed of light
%invert line profile and normalise
int=1-int;
int=int/max(int);
%function parameters
a1=locs;%input line centre (in km/s), this will be subtracted from the velocities to centre the line
vel=vel-a1;
a3=10;%input starting V sin(i) (in km/s, a range will be tested)
a2=a3/2;%input starting FWHM for pre-rotationally broadened line (in km/s, a range will be tested)
lamda0=5500;%input lambda0 (use centre wavelength of ccf if a ccf)
save('temp_lambda0','lamda0')
lb=[0 0];
ub=[200 1000];
%respfnwidth=input('input spectrograph response width parameter or leave empty to fit Thar lines (telwave and telint); -> ');
epsx=0.555;%input limb darkening coefficient
save('temp_epsx','epsx')
%convert FWHM to wavelength
x0=[a2,a3];
x0(1)=a2(1)*lamda0/c;
%fit line data with a broadened Gaussian
opts = optimset('Display','off');
%figure
[x,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@vsini_fittingfn,x0,vel,int,lb,ub,opts);
results=x;
x(1)=x(1)/lamda0*c;
results1=x;
width=x(1); % for use later
FWHMgauss=x(1);
clear x0 x
vsini=results1(2);
fittedint=vsini_fittingfn(results,vel);
fittedint=1-fittedint;
int=1-int;
% figure
% plot(vel,1-int,'.k')
% hold on
% plot(vel,1-fittedint,'r')
% hold off
% xlabel('Radial velocity (km/s)')
% ylabel('Intensity')
% legend('Observed line','Fitted line','Location','SouthEast')
% axis([-inf inf -0.05 1.1])
% disp(' ')
% fprintf('\nResults:\nFWHM for line and spectrograph response fn = %.2f km/s\nV sin(i) = %.2f km/s\n\n',results1(1),results1(2));
