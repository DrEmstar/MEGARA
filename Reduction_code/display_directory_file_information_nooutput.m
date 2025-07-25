function info = display_directory_file_information_nooutput()
%This code queries the J*.fit files in a directory to dertermine their basic contents without a display

file_list=dir('J*.fit');
fid=fopen('file_information.txt','w');
for s=1:numel(file_list)
    nam=file_list(s).name;nam(end-3:end)=[];
    eval(['header = get_header(''' file_list(s).name ''');'])
    temp = get_required_headers_from_header(header);
    temp=check_hercules_headers(temp);
    eval(['info.list{s}=''' nam ''';'])
    eval(['info.' nam '=temp;'])
    fprintf(fid,'%s  %7s   %.4f    %d    %7.2f     %s\n',file_list(s).name,temp.HERCEXPT,temp.MJD_OBS,temp.HERCFIB,temp.EXPTIME,temp.OBJECT);
end
fclose(fid);
if isempty(file_list)
    info.list=[];
end

