function [backfit,yy,xx,zz] = fit_background(summed_flat,allorders)
%select points for background fit

test=95-allorders.numords;
for s=1:allorders.numords
    if s<=90-test
        m=151-s-test;
        Q(s,:)=[s m];
    else
        m=150-s-test; %due to missing order (m=60)
        Q(s,:)=[s m];
    end
end
m=Q(:,2);

xx=1:size(summed_flat,2);
st=160;
numords=allorders.numords;
cnt=0;

for D=2:numords
    cnt=cnt+1;
    eval(['minpoint=max( [min(allorders.points' num2str(D-1) ') min(allorders.points' num2str(D) ')] );'])
    eval(['maxpoint=min( [max(allorders.points' num2str(D-1) ') max(allorders.points' num2str(D) ')] );'])

    x=minpoint:st:maxpoint;
    eval(['S=find(allorders.points' num2str(D-1) '==min(x));'])
    eval(['SSS=find(allorders.points' num2str(D) '==min(x));'])
    s=S:st:S+st*(numel(x)-1);
    sss=SSS:st:SSS+st*(numel(x)-1);
    
   %Note from Emily- I am not sure of the function of this code- to check one day. 
    if Q(D,2)<62
        if Q(D,2)>58.5 && Q(D,2)<59.5
            eval(['ty=diff( cat(1,allorders.ofit' num2str(D-1) '(s),allorders.ofit' num2str(D) '(sss)) );'])
            eval(['y=round(allorders.ofit' num2str(D-1) '(s)+ty/4);'])
        else
            eval(['y=round(mean(cat(1,allorders.ofit' num2str(D-1) '(s),allorders.ofit' num2str(D) '(sss))));'])
        end
        pos=y(find_nearest(x,size(summed_flat,2)/2));
        [pp,S,mu] = polyfit(x,y,3);
        yy(cnt,:)=polyval(pp,xx,[],mu);
        lastyy=yy(cnt,:);
        posgood=lastyy( round(numel(lastyy)/2) );
        %startfrom goodind
        goodind=Q(Q(:,2)==62,1);
        yy(cnt,:)=yy(cnt-1,:)+( pos-posgood );
    else
        eval(['y=round(mean(cat(1,allorders.ofit' num2str(D-1) '(s),allorders.ofit' num2str(D) '(sss))));'])
        [pp,S,mu] = polyfit(x,y,3);
        yy(cnt,:)=polyval(pp,xx,[],mu);
    end
end

T=medfilt2(summed_flat,[5 5]);
%extend points to either side of end orders
y1=yy(1,:)-( diff(yy(1:2,:)) );
yend=yy(end,:)+( diff(yy(end-1:end,:)) );
yy=cat(1,y1,yy,yend);

yy=round(yy);
zz=zeros(size(yy));
for n=1:size(yy,1)
    for m=1:size(yy,2)
        if yy(n,m)>=1 && yy(n,m)<=size(summed_flat,1);
            zz(n,m)=T(yy(n,m),xx(m));
        end
    end
end

xx=repmat(xx,[size(yy,1) 1]);
[Y,X] = meshgrid(1:size(summed_flat,1),1:size(summed_flat,2));

%interpolating background points
warning('off','MATLAB:griddata:DuplicateDataPoints');

%remove nans
yy(isnan(yy))=1;

backfit = griddata(yy(:),xx(:),zz(:),Y,X,'cubic');

B=isnan(backfit);
backfit(B)=0;
backfit=backfit';