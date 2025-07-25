function remove_files(badfiles)

warning off
mkdir('removed_files');
warning on

movepath=strcat(cd,'/removed_files');

%move fit files to  folder
for n=1:numel(badfiles)
    file = char(badfiles(n));
    try
        movefile(file, movepath);
    catch
        disp(strcat(file, ' could not be moved.'))
        continue
    end
end
files=dir('g*.mat');
for n=1:numel(files)
    delete(strcat(cd, '/', files(n).name))
end
%remove mat files
files2=dir('bad*J*.mat');
for n=1:numel(files2)
    delete(strcat(cd, '/', files2(n).name))
end
end