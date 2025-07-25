function [raw_data_location,reduced_data_location]=define_directory_locations
%setting up where all the directories are
raw_data_directory=which('pointer_raw_data_directory.m');
if isempty(raw_data_directory)==1
    disp('No raw data directory found. Move the pointer file to the data directory and check the directory is added to the set path')
    return
end
raw_data_location=raw_data_directory(1:end-28);

reduced_data_directory=which('pointer_reduced_data_directory.m');
if isempty(reduced_data_directory)==1
    disp('No reduced data directory found. Move the pointer file to the data directory and check the directory is added to the set path')
    return
end
reduced_data_location=reduced_data_directory(1:end-32);