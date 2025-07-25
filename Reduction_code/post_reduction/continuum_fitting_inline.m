function [result,thecfit] = continuum_fitting_inline(wave,int,synthint)

cfit=[];
[thecfit,result] = update_plot(wave,int,cfit,synthint);

while 1
    disp('1 - add point')
    disp('2 - delete point')
    disp('3 - move point')
    disp('4 - exit')
    disp(' ')
    entry = input('-> ');
    if isempty(entry)
        entry=0;
    end
    switch entry
        case 1
            [result,thecfit,cfit]=add_point(wave,int,cfit,synthint);
        case 2
            [result,thecfit,cfit]=delete_point(wave,int,cfit,synthint);
        case 3
            [result,thecfit,cfit]=move_point(wave,int,cfit,synthint);
        case 4
            break
        otherwise
            disp('not one of the available options ...')
    end
    disp(' ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,thecfit,cfit]=add_point(wave,int,cfit,synthint)
zoom off
pan off
figure(100)
disp('waiting for mouse button press within axes')
x=ginput(1);
cfit=sortrows(cat(1,cfit,x),1);
try
    [thecfit,result]=update_plot(wave,int,cfit,synthint);
catch
    disp('error in update plot, possibly two points at same place, try again!')
    thecfit=[];
    result=[];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,thecfit,cfit]=delete_point(wave,int,cfit,synthint)
if isempty(cfit)
    disp('no points to delete!')
    thecfit=[];
    result=[];
    return
end
zoom off
pan off
figure(100)
disp('waiting for mouse button press within axes')
xy_values=ginput(1);
xy_nearest = findnearest2D(cfit,xy_values);
cfit(xy_nearest,:)=[];
try
    [thecfit,result]=update_plot(wave,int,cfit,synthint);
catch
    disp('error in update plot, possibly two points at same place, try again!')
    thecfit=[];
    result=[];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,thecfit,cfit]=move_point(wave,int,cfit,synthint)
if isempty(cfit)
    disp('no points to move!')
    thecfit=[];
    result=[];
    return
end
zoom off
pan off
disp('waiting for mouse button press within axes')
figure(100)
while 1
    w = waitforbuttonpress;
    if w == 0
        pos=get(gca,'CurrentPoint');
        pos=pos(1,1:2);
        xy_nearest = findnearest2D(cfit,pos);
        plot(cfit(xy_nearest,1),cfit(xy_nearest,2),'ko')
        disp('circled point will move to position of next mouse button press ...')
        break
    end
end
disp('waiting for mouse button press within axes')
figure(100)
while 1
    w = waitforbuttonpress;
    if w == 0
        pos=get(gca,'CurrentPoint');
        pos=pos(1,1:2);
        cfit(xy_nearest,:)=pos;
        try
            [thecfit,result]=update_plot(wave,int,cfit,synthint);
        catch
            disp('error in update plot, possibly two points at same place, try again!')
            thecfit=[];
            result=[];
            return
        end
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thecfit,result] = update_plot(wave,int,cfit,synthint)
figure(100)
axes(1)=subplot(2,1,1,'replace');
hold on
title('Input and continuum fit')
ylabel('Intensity (unknown)')
xlabel('Wavelength')
plot(wave,int,'Parent',axes(1))
plot(wave,int./synthint,'g','Parent',axes(1))

if numel(cfit) > 0
    plot(cfit(:,1),cfit(:,2),'+r','Parent',axes(1))
end
if length(cfit) >= 7
    thecfit=spline(cfit(:,1),cfit(:,2),wave);
    result=int./thecfit;
    plot(wave,thecfit,'g','Parent',axes(1))
    axes(2)=subplot(2,1,2,'replace');
    hold on
    title('Resulting fitted data')
    xlabel('Wavelength (Angstroms)')
    ylabel('Normalised Intensity')
    axis([-inf inf -inf inf])
    linkaxes(axes,'x');
    plot(wave,synthint,'r','Parent',axes(2))
    plot(wave,result,'Parent',axes(2))
else
    thecfit=[];result=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xy_nearest = findnearest2D(vector2D,xy_values)
sz=size(vector2D);
if sz(2) > 2
    disp('problem with input vector for findnearest2D')
    return
end
dist=sqrt( (vector2D(:,1)-xy_values(1)).^2 + (vector2D(:,2)-xy_values(2)).^2 );
xy_nearest=find(dist==min(dist));
if numel(xy_nearest) > 1
    xy_nearest=xy_nearest(floor(numel(nearest)/2));
end
