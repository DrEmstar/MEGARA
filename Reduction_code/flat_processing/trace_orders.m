function [allorders, summed_flat] = trace_orders(summed_flat,minpixseparation)

sz=size(summed_flat);
x=1:sz(1);
tracedata=cat(2,x',mean(summed_flat(:,sz(2)/2-5:sz(2)/2+5)')');
tracedata(:,2)=tracedata(:,2)-min(tracedata(:,2));
S=running_min(tracedata(:,2),10);
S=mean_smoothing(medfilt1(S,21),10);
tracedata(:,2)=tracedata(:,2)-S;
tracedata(:,2)=tracedata(:,2)./mean_smoothing(running_max(tracedata(:,2),30),30);

X=find(tracedata(:,2)>0.25);
nn=find(diff(X) > minpixseparation);
p=tracedata(X(nn),1);p(1)=[];
q=tracedata(X(nn+1),1);q(end)=[];
orders=mean(cat(2,p,q)');

disp('getting starting centre orders ...')

%fit 1-D Gaussian at 20 pixels around original trace
width=0.4;
fitsize=15;
centre=zeros(numel(orders),1);
widths=centre;
sse=centre;
Fits=zeros(length(round(orders(1)-tracedata(1)+1-fitsize):round(orders(1)-tracedata(1)+1+fitsize)) ,numel(orders)*3);
for s=1:numel(orders)
    try
        fitdata=tracedata(round(orders(s)-tracedata(1)+1-fitsize):round(orders(s)-tracedata(1)+1+fitsize),:);
    catch
        fitdata=tracedata(round(orders(s)-tracedata(1)+1-fitsize):numel(tracedata(:,1)),:);
        Fits=Fits(1:numel(fitdata(:,1)),:);
    end
    Fits(:,s*3-2:s*3-1)=fitdata;
    [centre(s), widths(s), sse(s), Fits(:,s*3)] = fit_lin_gauss(fitdata(:,1), fitdata(:,2));
end

S=sse < 0.2 & abs(widths) < 10 & abs(widths) > 2.5;
centre=centre(S);
widths=widths(S);
sse=sse(S);

%now we have starting centres for orders
%run an order tracing program that finds the order starting from the centre

disp('tracing each order ...')

ordcnt=0;
for s=1:numel(centre)
    [opos, ocentres, owidths]=trace_single_order_from_centre(summed_flat,centre(s),widths(s));
    S=owidths < 2.5 | owidths > 9;
    opos(S)=[]; ocentres(S)=[]; owidths(S)=[];

    if numel(opos)>10
        ordcnt=ordcnt+1;
        
        points=min(opos):max(opos);
        [p,S,mu]=polyfit(opos,ocentres,5);
        [ofit,delta]=polyval(p,points,S,mu);
        
        eval(['allorders.width' num2str(ordcnt) '=median(owidths);'])
        eval(['allorders.ofit' num2str(ordcnt) '=ofit;'])
        eval(['allorders.points' num2str(ordcnt) '=points;'])
    else
    end
end

allorders.numords=ordcnt;
