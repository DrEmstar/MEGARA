function [newwave] = shift_velocity(wave, vel)

c=299792.458; %speed of light (km/s)
shift=log(1+vel/c); %compute shift in log space
newwave=exp(log(wave)-shift); %create new shifted axis in wavelength space
