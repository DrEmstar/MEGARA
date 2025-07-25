function headers = get_required_headers_from_header(header)
headers.NAXIS1=[];
headers.NAXIS2=[];
headers.HERCEXPT=[];    %Hercules Exposure type
headers.HERCFIB=[];     %Hercules Fibre 1/2/3
headers.DATE_OBS=[];    %Observation date OLD CCD SOFTWARE
headers.DATE=[];        %Observation date NEW CCD SOFTWARE 
headers.MJD_OBS=[];     %Modfd, corr'd, heliocentric, mid obs JD
headers.REC_STRT=[];    %Recorded UT time of observation start
headers.EXPTIME=[];     %Actual exposure time in seconds
headers.HERCFWMT=[];    %Hercules flux weighted mean exp time (mins) ->actually seconds I think ...
headers.HERCFTC=[];     %Hercules flux meter total counts
headers.PORT=[];        %number of ports -1
headers.RA=[];          %RA J2000 coords
headers.DEC=[];         %DEC J2000 coords
headers.OBSERVER=[];
%following are HERCULES processed headers (using HRSP: prepare_descriptors)
headers.RVCORR=[];      %
headers.JDCORR=[];      %
headers.JD_MID=[];      %
headers.OBJECT=[];      %


for m=1:numel(header)
    if ~isempty(findstr(header{m},'NAXIS1'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.NAXIS1=t;
    elseif ~isempty(findstr(header{m},'NAXIS2'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.NAXIS2=t;
    elseif ~isempty(findstr(header{m},'HERCEXPT'))
        t=header{m};t=t(10:31);
        t(double(t)==32)=[]; %trim off spaces
        t(double(t)==39)=[]; %trim off apostrophes
        if contains(t,'Thorium') || contains(t,'THORIUM')
            headers.HERCEXPT='Thorium';
        elseif contains(t,'Stellar') || contains(t,'STELLAR')
            headers.HERCEXPT='Stellar';
        elseif contains(t,'White') || contains(t,'WHITE') || contains(t,'Smooth') || contains(t,'SMOOTH')
            headers.HERCEXPT='White L';
        end
    elseif ~isempty(strfind(header{m},'HERCFIB'))
        t1=header{m};t=t1(12);
        if t==' '
           t=t1(28);
        end
        t=str2double(t);
        headers.HERCFIB=t;
    elseif ~isempty(strfind(header{m},'DATE-OBS'))    
        t=header{m};t=t(10:31);t=strcat(t);t=t(2:end);
        headers.DATE_OBS=t;
    elseif ~isempty(strfind(header{m},'DATE'))         % NEW CCD SOFTWARE
        t=header{m};t=t(10:31);t=strcat(t);t=t(2:end); % NEW CCD SOFTWARE
        headers.DATE=t;                                % NEW CCD SOFTWARE
    elseif ~isempty(strfind(header{m},'MJD-OBS'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.MJD_OBS=t;
    elseif ~isempty(strfind(header{m},'REC-STRT'))
        t=header{m};t=t(10:31);t=strcat(t);t=t(2:end);
        headers.REC_STRT=t;
    elseif ~isempty(strfind(header{m},'EXPTIME'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.EXPTIME=t;
    elseif ~isempty(strfind(header{m},'HERCFWMT'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.HERCFWMT=t;
    elseif ~isempty(strfind(header{m},'HERCFTC'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.HERCFTC=t;
%   elseif ~isempty(strfind(header{m},'PARAM59')) 
    elseif ~isempty(strfind(header{m},'PARAM58')) % for NEW CCD SOFTWARE   
        t=header{m};t=t(10:30);
        t=str2double(t);
        if t>10% for OLD CCD SOFTWARE
            n=m+1;
            t=header{n};t=t(10:30);
            t=str2double(t);
            headers.PORT=t+1;
        else 
            headers.PORT=t+1;
        end
    elseif ~isempty(strfind(header{m},'RVCORR'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.RVCORR=t;
    elseif ~isempty(strfind(header{m},'JDCORR'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.JDCORR=t;
    elseif ~isempty(strfind(header{m},'JD-MID'))
        t=header{m};t=t(10:31);
        t=str2double(t);
        headers.JD_MID=t;
    elseif ~isempty(strfind(header{m},'OBSERVER'))
        t=header{m};t=t(10:31);
        headers.OBSERVER=t;
    elseif ~isempty(strfind(header{m},'COMMENT  RA'))
        t=header{m};
        try t=t(15:24);headers.RA_2000=t;
        catch headers.RA_2000=[];
        end
    elseif ~isempty(strfind(header{m},'COMMENT  DEC'))
        t=header{m};
        try t=t(16:24);headers.DEC_2000=t;
        catch headers.DEC_2000=[];
        end
    elseif ~isempty(strfind(header{m},'OBJECT'))
        try
            try
                %trim off extras
                t=header{m}; t=t(10:30); 
            catch
                %trim off extras
                t=header{m}; t=t(10:end); 
            end
            %trim off spaces
            t(double(t)==32)=[]; 
            %trim off apostrophes
            t(double(t)==39)=[]; 
            headers.OBJECT=t;
        catch
            disp('no OBJECT header found!')
        end
    end
end

%check for important missing headers!
if isempty(headers.MJD_OBS)
    if isempty(headers.JD_MID)
        if isempty(headers.DATE_OBS) && isempty(headers.EXPTIME) && isempty(headers.REC_STRT)
            disp('NOTE: there is insufficient header information to determine JD of observation!!!!!!!')
        end
    end
end

