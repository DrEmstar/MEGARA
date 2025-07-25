function data = raw_hercules_fitsread(filename)

header = get_header(filename);
headers = get_required_headers_from_header(header);

fid=fopen(filename,'r','ieee-be');
data=fread(fid,'int16');
fclose(fid);
clear fid ans

n=658; %the offset for reading data (including making bad start and end rows all appear at start)

%make reading offset disappear
data=data-min(min(data)); 
data(1:numel(data)-headers.NAXIS1*headers.NAXIS2-n)=[];
data(end-n+1:end)=[];
try
    data=reshape(data,[headers.NAXIS1 headers.NAXIS2 1]);
catch
    disp('file could not be read')
    badfile={filename};
    remove_files(badfile);
    return
end
%trimming poor areas of the image
data(1:34,:)=[];   %top and bottom bad parts (read in to be all at the start)
data(:,end-29:end)=[];  %edges
data(:,1:30)=[];

%rotate image so that it is blue to red bottom to top and left to right
data=fliplr(data);
