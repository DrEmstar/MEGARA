function [pos, centres, widths] = trace_single_order_from_centre(summed_flat,cen,wid)

step=10;
failtimes=20;

failcheck(1:failtimes)=1;
sz=size(summed_flat);
%trace from centre to left in 'step' pixel steps
imglrcen=round(sz(2)/2); %centre pixel
lp=imglrcen:-step:0;
if lp(end)==0
    lp(end)=1;
elseif lp(end)~=1
    lp(end+1)=1;
end
centresl=zeros(numel(lp),1);
widthsl=centresl;

cnt=0;
failcnt=0;
for s=lp
    cnt=cnt+1;
    if s==imglrcen
        newcen=cen;
    end
    if round(newcen-3*wid) < 1 || round(newcen+3*wid) > sz(1) %then outside image so exit
        break
    end
    xdata=round(newcen-3*wid:newcen+3*wid);
    ydata=mean(summed_flat(round(newcen-3*wid:newcen+3*wid),s:s+3)')';
    ydata=ydata(:)-min(ydata);
    [centresl(cnt), widthsl(cnt)] = find_max_and_fwhm(xdata(:), ydata);
    
    if abs(centresl(cnt)-newcen)>5 || widthsl(cnt)<2.5 || widthsl(cnt)>11
        %if fail failtimes consecutive times then exit
        failcnt=failcnt+1;
        centresl(cnt)=0; widthsl(cnt)=0;
        fail(failcnt)=cnt;
        e=diff(fail);
        if numel(e)>=failtimes
            if e(end-(failtimes-1):end)==failcheck
                break
            end
        end
    end
    newcen=centresl(cnt);
end
SS=centresl==0;
lpos=lp;
lpos=lpos(:);
lpos(SS)=[]; centresl(SS)=[]; widthsl(SS)=[];
%now go centre to right
rp=imglrcen+step:step:sz(2);
if rp(end)~=sz(2)
    rp(end+1)=sz(2);
end
centresr=zeros(numel(rp),1);
widthsr=centresr;
cnt=0;
failcnt=0;
fail=[];
for s=rp
    cnt=cnt+1;
    if s==imglrcen+step
        newcen=cen;
    end
    if round(newcen-3*wid) < 1 || round(newcen+3*wid) > sz(1) %then outside image!
        continue
    end
    xdata=round(newcen-3*wid:newcen+3*wid);
    ydata=mean(summed_flat(round(newcen-3*wid:newcen+3*wid),s-3:s)')'; ydata=ydata(:)-min(ydata);
    [centresr(cnt), widthsr(cnt)] = find_max_and_fwhm(xdata(:), ydata);
    if abs(centresr(cnt)-newcen)>5 || widthsr(cnt)<2.5 || widthsr(cnt)>11
        %if fail failtimes consecutive times then exit
        failcnt=failcnt+1;
        centresr(cnt)=0; widthsr(cnt)=0;
        fail(failcnt)=cnt;
        e=diff(fail);
        if numel(e)>=failtimes
            if e(end-(failtimes-1):end)==failcheck
                break
            end
        end
    end
    newcen=centresr(cnt);
end
SS=centresr==0;
rpos=rp;rpos=rpos(:);
rpos(SS)=[]; centresr(SS)=[]; widthsr(SS)=[];
%arrange the fits
centresl=flipud(centresl); widthsl=flipud(widthsl); lpos=flipud(lpos);
pos=cat(1,lpos,rpos);
centres=cat(1,centresl,centresr);
widths=cat(1,widthsl,widthsr);

