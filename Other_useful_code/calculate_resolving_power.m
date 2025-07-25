function [resolving_power,wave_loc, number_lines]=calculate_resolving_power(thar_image,flat_info,blue_data_chop,thar_info)

%Calculation of resolving power of spectral observations from ThAr frames

load(flat_info)
thar_info_data=load(thar_info);
img=hercules_fitsread(thar_image);

%get the wavelength solution
head=get_header(thar_image)';
header=get_required_headers_from_header(head);
header=check_hercules_headers(header);
jd=header.MJD_OBS;

thar_index=find_nearest(thar_info_data.jd_ths,jd);
fwave=eval(['thar_info_data.fwave.th' num2str(thar_index)]);

img(1:blue_data_chop,:)=[];

extracted_thar= extract_all_orders_no_background(allorders,img);

for s=1:flatdata.numords
    eval(['thardata.order_' num2str(s) '=sum(extracted_thar.order_' num2str(s) '.data(:,flatdata.order_' num2str(s) '.extraction_inds),2);'])
end
thardata=redo_order_numbering_simple(thardata,allorders);

figure
hold on
for m=1:numel(thar_info_data.m)
    ord=thar_info_data.m(m);
    %locate peaks of thars
    eval(['[pks, locs, w]=findpeaks(thardata.order_' num2str(ord) '-min(thardata.order_' num2str(ord) '),fwave.order_' num2str(ord) ',''Annotate'',''extents'',''MinPeakProminence'',20000,''WidthReference'',''halfheight'');']);  
    %cut lines and fit shoulders for a better fit.
    for n=1:numel(pks)
        pkextent=w(n);
        pk_left=locs(n)-pkextent;
        pk_right=locs(n)+pkextent;
        eval(['cut_line_fwave=fwave.order_' num2str(ord) '(find_nearest(fwave.order_' num2str(ord) ',pk_left):find_nearest(fwave.order_' num2str(ord) ',pk_right));']);
        eval(['cut_line_int=thardata.order_' num2str(ord) '(find_nearest(fwave.order_' num2str(ord) ',pk_left):find_nearest(fwave.order_' num2str(ord) ',pk_right));']);
        cut_line_int=cut_line_int-min(cut_line_int);
        [cut_pk(n), cut_locs(n), cut_w(n)]=findpeaks(cut_line_int,cut_line_fwave,'Annotate','extents','NPeaks',1,'WidthReference','halfheight');
    end
    res_powers=locs./w;
    keep_ind=find(res_powers>30000);
    res_powers=res_powers(keep_ind);
    eval(['resolving_power.order_' num2str(ord) '=sum(cut_pk(keep_ind).*res_powers'')./sum(cut_pk(keep_ind));']);
    eval(['number_lines.order_' num2str(ord) '=numel(keep_ind);']);
    eval(['wave_loc.order_' num2str(ord) '= sum(cut_pk(keep_ind).*cut_locs(keep_ind))./sum(cut_pk(keep_ind));']);
    eval(['plot(thar_info_data.m(m),resolving_power.order_' num2str(ord) ',''x'')'])
    keep_ind=[];
    res_powers=[];
    cut_pk=[];
    cut_locs=[];
    cut_w=[];
end
xlabel('Order Number')
ylabel('Resolving Power')
