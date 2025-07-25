function [type,star] = identify_the_file(files,s)
eval(['type=files.' files.list{s} '.HERCEXPT;']); % either 'Thorium','White L' or 'Stellar'.
eval(['star=files.' files.list{s} '.OBJECT;']); %'returns name of star for that file
end