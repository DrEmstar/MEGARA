function has_run=badfile_check_general
%see if this has already run- skip if not
has_run=0;
badfiles_general=[];
readout1=[];
readout2=[];
test=dir('badfiles_general.mat');
if numel(test)>0
    %disp('General file check already completed in this directory')
    has_run=1;
    return
end

%remove duplicate files sometimes generated
badlist=dir('J*@*.fit');
if isempty(badlist)==0
    for k = 1 : length(badlist)
        delete(badlist(k).name);
    end
end

filesize = dir('J*.fit');
%this checks file size. anything <33.8mb will be bad!!!!
for s=1:numel(filesize)
    if filesize(s).bytes <33842880
        badfiles_general=  cat(1,badfiles_general, [{strcat(cd,'/', filesize(s).name)}]);
        readout1 = cat(1,readout1, [strcat(filesize(s).name, ' is an incomplete file.')]);
        continue
    end
    if filesize(s).bytes >34000000
        badfiles_general=  cat(1,badfiles_general, [{strcat(cd,'/', filesize(s).name)}]);
        readout2 = cat(1,readout2, [strcat(filesize(s).name, ' is a large file.')]);
        continue
    end
end
remove_files(badfiles_general)
info = display_directory_file_information_nooutput();
save('file_info','info')
if has_run==0
    save('badfiles_general', 'badfiles_general')
    has_run=1;
end
