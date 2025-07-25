function header=check_hercules_headers(head)
header=head;
try
    if isempty(head.MJD_OBS)
        if ~isempty(head.DATE_OBS) & ~isempty(head.REC_STRT) & ~isempty(head.EXPTIME)
            date=head.DATE_OBS;
            year=str2num(date(2:5));month=str2num(date(7:8));day=str2num(date(10:11));
            time=head.REC_STRT;time(1)=[];time(end)=[];time(end+1)='0';
            JDstart = convert_date_to_JD(year,month,day,time);
            %flux weighted exposure time added
            expt=head.HERCFWMT/86400;
            JDmid=JDstart+expt;
            %to keep inline with normally written fits headers
            JDmid=JDmid-2.4e6;
            header.MJD_OBS=JDmid;
        else
            disp('not enough information to reconstruct MJD time from DATE_OBS REC_STRT and EXPTIME, exiting ...')
            return
        end
    end
catch
    return
end
