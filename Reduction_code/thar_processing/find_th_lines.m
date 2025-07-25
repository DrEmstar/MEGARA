function [thfitinfo,numbad_thar_flag] = find_th_lines(thdata,order,ul,uc,ur,plot_q)

cnt=0;
cnt2=0;
d=numel(uc);
for t=1:d
    try
        eval(['xdata=thdata.order_' num2str(order(t)) '.xax;'])
        eval(['ydata=thdata.order_' num2str(order(t)) '.summed_data;'])
        xdata=xdata(:);
        ydata=ydata(:);       
    catch
        continue
    end
    if isempty(xdata)
        continue
    end
    if ul(t)>=min(xdata) && ur(t)<=max(xdata)
        cnt=cnt+1;
        S=xdata>=round(ul(t)) & xdata<=round(ur(t));
        xdata=xdata(S);
        ydata=ydata(S);
    else
        continue
    end
    if isempty(ydata)
        continue
    end
    
    ydata=ydata-min(ydata);
    if ~mean(ydata)==0
        ydata=ydata/mean(ydata);
    end
    
    lb = [0 ul(t) 1 -0.5 -1];
    startpoint = [1 uc(t) 4 0 0];
    ub = [50 ur(t) 15 1.0 1];
    
    [x,resnorm,residual,exitflag,fitted] = fit_th_line(xdata,ydata,lb,ub,startpoint);
    
    if isempty(plot_q) || numel(plot_q)~=1
        plot_q=0;
    end
    if plot_q==1
        fprintf('resnorm=%.3f\n',resnorm)
        fprintf('hght=%.2f posn=%.1f FWHM=%.2f pos_ud=%.4f slope=%.4f\n',x(1),x(2),x(3),x(4),x(5))
        disp(' ')    
        plot(xdata,ydata)
        title(['line number ' num2str(t)])
        hold on
        plot(xdata,fitted,'r')
        hold off
        tes=input('enter x to stop displaying ThAr fits -> ','s');
        if strcmp(tes,'x')
            plot_q=0;
        end
    end
    
    if resnorm<=2 && x(1)>0.8 && x(3)>=1.5 && x(3)<=9
        cnt2=cnt2+1;
        thfitinfo(cnt2,:)=cat(2,t,resnorm,x);
    end
end
if length(thfitinfo)<800
    numbad_thar_flag=1;
else 
    numbad_thar_flag=0;
end

fprintf('%.0f lines found\n',length(thfitinfo))


