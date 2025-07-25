function [wavelengths, ele, elem, widths] = read_in_output12

findfile=dir('output.12');
if numel(findfile) ~= 1
    disp('there must be a copy of the file output.12 in the current directory')
    return
end
fid = fopen('output.12');
C = textscan(fid, '%n %n %n %s %s %n %n %n %n %*s %n %n %n');
fclose(fid);
wavelengths=C{1,3};
ele=C{1,4};
ion=C{1,5};
widths=C{1,9};
elem=strcat(ele,ion);
d=find(widths == 0);
ele(d)=[];
elem(d)=[];
wavelengths(d)=[];
widths(d)=[];
