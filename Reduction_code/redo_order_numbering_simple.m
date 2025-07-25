function [Data m]=redo_order_numbering_simple(data,allorders)

%find the mean y-pixel position of each of the orders
for s=1:allorders.numords
    eval(['meanordypos(s) = mean(allorders.ofit' num2str(s) ');'])
end
%find where the found orders skips one, that skipped order is order 60
ind=find( diff(meanordypos(70:end)) == max(diff(meanordypos(70:end))) ); 
ind=ind+70; %the order (in initial numbering) that indicates order 59
orig_nums=1:allorders.numords;
m(1:ind-1)= 61+(ind-1)-1:-1:61;
m(ind:allorders.numords)=59:-1:59-(allorders.numords-ind);
for s=1:allorders.numords
    eval(['Data.order_' num2str(m(s)) '=data.order_' num2str(s) ';'])
end
