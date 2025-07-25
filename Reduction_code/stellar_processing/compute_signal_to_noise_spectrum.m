function signal_to_noise=compute_signal_to_noise_spectrum(stellar_order,exp_time)

load photon_ADU_conversion.mat

%dark noise
dark=find_dark_value(exp_time);
load readout_noise_value.mat
dark=dark-bias;

%shotnoise- extracted order 100 from the stellar image
num_pix=numel(stellar_order.data);
img_sum=(sum(stellar_order.data,2)/size(stellar_order.data,2));

%scaling backgrounds
bias=bias/(4096*4036);
sigma_bias=sqrt(bias);
dark=dark/(4096*4036);
sigma_dark=sqrt(dark);

sigma_shot=sqrt(img_sum-bias-dark);

%backround total
B=sqrt(sigma_shot.^2+sigma_dark^2+sigma_bias^2)*mu_0_avg;

%light total
N=img_sum*mu_0_avg;
L=(N-B);

signal_to_noise_spectrum=L./(sqrt(L+2*B));

%converting to a single value using the section of order 100 where there is no overlap
signal_to_noise_spectrum=signal_to_noise_spectrum(1700:2219);
signal_to_noise=median(signal_to_noise_spectrum);
signal_to_noise=round(signal_to_noise,0);

