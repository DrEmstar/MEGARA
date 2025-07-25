function [synthwave, synthint] = read_in_output11

findfile=dir('output.11');
if numel(findfile) ~= 1
    disp('there must be a copy of the file output.11 in the current directory')
    return
end

fid = fopen('output.11');
C = textscan(fid, '%n %n');
fclose(fid);
synthwave=C{1};
synthint=C{2};
