function [jd_ths, m, fwave, numbad_thar,save_flag] = thar_processing(info, allorders, Bias, flatdata, blue_data_chop_value,numbad_thar)

save_flag=1;
Thar_info_name=['ThAr_info_blue_' num2str(blue_data_chop_value) '.mat'];
test=dir(Thar_info_name);
if numel(test)>0
    load(Thar_info_name)
    disp('ThAr_info loaded from previous run')
    save_flag=0;
    %check successful ThArs in fwave
    cntt=0;
    while 1
        cntt=cntt+1;
        try
            eval(['tt=fwave.th' num2str(cntt) ';'])
        catch
            break
        end
    end
    cnttt=0;
    %check number of ThArs in the directory according to info
    for xx=1:numel(info.list)
        eval(['temp=info.' info.list{xx} '.HERCEXPT;'])
        if strcmp(temp,'Thorium') || strcmp(temp,'THORIUM')
            cnttt=cnttt+1;
        end
    end
    cnt=0;
    for s=1:numel(info.list)
        eval(['temp=info.' info.list{s} '.HERCEXPT;'])
        if strcmp(temp,'Thorium') || strcmp(temp,'THORIUM')
            cnt=cnt+1;
            if cnt >= cntt
                nam=[info.list{s} '.fit'];
                fprintf('processing ThAr spectrum file: %s\n',nam)
                eval(['[jd_ths(cnt), m, fwave.th' num2str(cnt) ']=process_thar_frames(nam, info, allorders, Bias, flatdata, blue_data_chop_value);'])
            end
        end
    end
    return
end

parfor s=1:numel(info.list)%parfor
    [temp1,temp2,temp3] = make_temp_variable_thars(info,s);
    if strcmp(temp1,'Thorium') || strcmp(temp1,'THORIUM')
        %if temp3<12 && temp3>3
            nam=[info.list{s} '.fit'];
            fprintf('processing ThAr spectrum file: %s\n',nam)
            [jd_ths_1, m_temp, fwave_temp,numbad_thar_flag]=process_thar_frames(nam, info, allorders, Bias, flatdata, blue_data_chop_value);
            jd_ths_t(s)=jd_ths_1;
            numbad_thar_list(s)=numbad_thar_flag;
            fwave_t{s}= fwave_temp;
            m_t{s}=m_temp;
        %else
        %end
    end
end

numbad_thar=numel(find(numbad_thar_list));

test_jd_ths_t=exist('jd_ths_t','var');
if test_jd_ths_t==1
    indie=find(jd_ths_t);
    jd_ths=jd_ths_t(indie);
    for n=1:numel(indie)
        eval(['fwave.th' num2str(n) '=fwave_t{indie(n)};'])
        if n==1
            m=m_t{indie(n)};
        end
    end
end