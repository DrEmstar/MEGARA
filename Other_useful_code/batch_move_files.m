%filter out files with certain parameters. This example filters all flats with exposure times greater than 3 seconds.
clear all
clc

load file_info.mat
files=info.list;
for n=1:numel(files)
    file=info.list{n};
    exp_type=eval(['info.' file '.HERCEXPT']);
    exp_time=eval(['info.' file '.EXPTIME']);
    if strcmp(exp_type,'Thorium') && exp_time>6
        eval(['movefile(''' file '.fit'',''removed_files'')'])
    end
end


