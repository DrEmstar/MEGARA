function stardata_c = remove_cosmics(stardata,plotting)

try
    numords=stardata.numords;
catch
    disp('Can''t determine number of orders -> skipping')
    stardata_c=stardata;
    return
end
cosmic_count=0;

for ord=1:numords
    try
        eval(['yax=stardata.order_' num2str(ord) '.yax;'])
        eval(['odata=stardata.order_' num2str(ord) '.data;'])
        eval(['xax=stardata.order_' num2str(ord) '.xax;'])
        eval(['ypos=stardata.order_' num2str(ord) '.ypositions;'])
        xax=repmat(xax(:),1,size(odata,2));
        odatab=odata;
        O=medfilt2(odata,[7 1]);
        ODATA=odata./mean(mean(odata));
        
        RES = ODATA - medfilt2(ODATA,[7 1]);
        C=RES>0.5;
        [x,y]=find(C);
        cosmic_count=cosmic_count+numel(x);
        %fix up odatab
        odatab(C) = O(C);
        
        eval(['stardata.order_' num2str(ord) '.data=odatab;'])
        eval(['stardata.order_' num2str(ord) '.data_pre_c=odata;'])
        eval(['stardata.order_' num2str(ord) '.cosmics=C;'])

        if plotting==1
            % commands to produce test spectrum
            figure(100)
            mesh(odatab)
            title(['cosmic ray fixed, order' num2str(ord)])
            figure(101)
            mesh(odata)
            title(['before cosmic ray fix, order' num2str(ord)])
            pause
        end
    catch
        fprintf('skipping order %.0f\n',ord)
    end
end

stardata_c=stardata;
