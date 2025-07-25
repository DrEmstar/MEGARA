function header = get_header(filename)

fid=fopen(filename,'r','ieee-be');
A=fread(fid,3*2880,'char');  
warning off
for x=1:36*3
    header{x}=strcat(char(A((x-1)*80+1:(x-1)*80+80))');
end
warning on
status=fclose(fid);

