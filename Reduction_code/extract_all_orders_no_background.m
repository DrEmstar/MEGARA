function extracted = extract_all_orders_no_background(allorders,img)

%get order numbers
try
    orders=1:allorders.numords;
catch
    disp('order numbers not defined in allorders.numorders ... exiting')
    extracted =[];
    return
end

%sequentially extract each order
cnt=0;
for ordnum=orders
    cnt=cnt+1;
    %get the things needed for each order - if they're not all there then skip the order
    try
        eval(['points=allorders.points' num2str(ordnum) ';'])
        eval(['ofit=allorders.ofit' num2str(ordnum) ';'])
        eval(['width=allorders.width' num2str(ordnum) ';'])
    catch
        fprintf('\nSkipping order %d as fitting information is not complete \n\n',ordnum)
        continue
    end
    
    %determine if WIDTHFACTOR HAS BEEN DETERMINED OR IF A STANDARD ONE SHOULD BE USED
    try
        widthfactor=allorders.WIDTHFACTOR;
    catch
        widthfactor=1.3; %standard value
    end
    
    try
        %extract the order
        eval(['order' num2str(ordnum) '= extract_order_no_background(points,ofit,width,img,widthfactor);'])
        %add extracted order to the 'extracted' output structure
        eval(['extracted.order_' num2str(ordnum) '=order' num2str(ordnum) ';'])
    catch
        fprintf('order %d has a problem and will not be extracted\n',ordnum)
    end
    
end

