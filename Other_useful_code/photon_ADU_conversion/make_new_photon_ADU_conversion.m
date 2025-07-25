%make_new_photon_ADU_conversion

%flats from tutorial data
%files={'J8480001.fit','J8480002.fit','J8480003.fit','J8480004.fit','J8480005.fit','J8480006.fit','J8480007.fit','J8480008.fit','J8480009.fit','J8480010.fit','J8479001.fit','J8479002.fit','J8479004.fit','J8479005.fit','J8479006.fit','J8479007.fit','J8479009.fit','J8479010.fit','J8477001.fit','J8477002.fit','J8477003.fit','J8477004.fit','J8477005.fit','J8477006.fit','J8477007.fit','J8477008.fit','J8477009.fit','J8477010.fit','J8476001.fit','J8476002.fit','J8476003.fit','J8476004.fit','J8476005.fit','J8476006.fit','J8476007.fit','J8476008.fit','J8476009.fit','J8476010.fit','J8464001.fit','J8464002.fit','J8464003.fit','J8464004.fit','J8464005.fit','J8464006.fit','J8464007.fit','J8464008.fit','J8464009.fit','J8464010.fit','J8478001.fit','J8478002.fit','J8478003.fit','J8478004.fit','J8478005.fit','J8478006.fit','J8467001.fit','J8467002.fit','J8467003.fit','J8466001.fit','J8466002.fit','J8466003.fit'};

files={'J7731002.fit','J7731003.fit','J7731004.fit','J7731005.fit','J7731007.fit'};

info = display_directory_file_information_nooutput();
Bias=0;
blue_data_chop_value=500;

%sum the flats first to get the orders defined
summed_flat = merge_flats(info,Bias);
summed_flat = blue_data_chop(summed_flat,blue_data_chop_value);
minpixseparation=10;
[allorders, summed_flat] = trace_orders(summed_flat,minpixseparation);

%get signal vector for order 100 for all flats
for n=1:numel(files)
    flat = raw_hercules_fitsread([files{n}]);
    flat = blue_data_chop(flat,blue_data_chop_value);
    extracted_flat = extract_all_orders_no_background(allorders,flat);
    flat_data=redo_order_numbering_simple(extracted_flat,allorders);
    flat_order=flat_data.order_100;
    num_pix=numel(flat_order.data);
    img_sum(n,:)=(sum(flat_order.data,2)/size(flat_order.data,2));
end

mean_flats=mean(img_sum,1);

%get noise vector for order 100 for all flats

for n=2:numel(files)
    f_scatter=img_sum(1,:)./img_sum(n,:);
    f_scatter(isinf(f_scatter)) = 1;
    rms_f_scatter(n-1,:)=rms(f_scatter,1);
end
mean_scatter=mean(rms_f_scatter);

mean_noise=mean_scatter.*mean_flats-mean_flats;

signal_to_noise=mean_flats./mean_noise;

load readout_noise_value.mat

bias_single_pix=bias/(4096*4036);

mu_0=(signal_to_noise.^2.*mean_flats)./(mean_flats-(bias_single_pix));
sigma=std(mu_0);
mu_0(mu_0>5*sigma)=median(mu_0);

mu_0_avg=median(mu_0(500:3500));

save('photon_ADU_conversion.mat','mu_0_avg')