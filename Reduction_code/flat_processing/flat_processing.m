function [allorders, flatdata, backdata, summed_flat, backfit, save_flag]=flat_processing(info,Bias,blue_data_chop_value,backdata_medfilt,plotting)

%Flat image processing
save_flag=1;
flatdataname=eval(['''flat_info_blue_' num2str(blue_data_chop_value) '.mat''']);
test=dir(flatdataname);
if numel(test)>0
    load(flatdataname)
    disp('flat_info loaded from previous run ...')
    save_flag=0;
    return
else
    summed_flat = merge_flats(info,Bias);
    summed_flat = blue_data_chop(summed_flat,blue_data_chop_value);
    minpixseparation=10;
    [allorders, summed_flat] = trace_orders(summed_flat,minpixseparation);

    disp('fitting the background ...')
    backfit= fit_background(summed_flat,allorders);
    summed_flat_b = summed_flat - backfit;
    disp('extracting the data ...')
    flatdata = extract_all_orders_no_background(allorders,summed_flat_b);
    disp('extracting the background data ...')
    backdata = extract_all_orders_no_background(allorders,backfit);
    
    %cosmic ray filtering, just in case
    flatdata.numords=allorders.numords;
    flatdata = remove_cosmics(flatdata,plotting);
    
    disp('finalising flat-field data')
    try
        for ord=1:allorders.numords
            eval(['odata=flatdata.order_' num2str(ord) '.data;'])
            ind=mean(odata) > max(mean(odata))/100; %mean contribution of that pixel line must be greater than 1% of maximum
            eval(['flatdata.order_' num2str(ord) '.extraction_inds=ind;'])
            eval(['flatdata.order_' num2str(ord) '.sumdata=sum(odata(:,ind),2)'';'])
            
            eval(['backdata.order_' num2str(ord) '.summed_data=sum(backdata.order_' num2str(ord) '.data(:,ind),2);'])
            eval(['backdata.order_' num2str(ord) '.summed_smoothed_data=mean_smoothing(backdata.order_' num2str(ord) '.summed_data,backdata_medfilt);'])
            
            eval(['flatdata.order_' num2str(ord) '.sumdata_b=flatdata.order_' num2str(ord) '.sumdata(:)-backdata.order_' num2str(ord) '.summed_smoothed_data(:);'])
        end
    catch
        disp('errors with finalising flat field data')
    end
end

disp('Flat-field images processing complete')
