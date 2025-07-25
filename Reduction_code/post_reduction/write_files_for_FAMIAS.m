function write_files_for_FAMIAS(jd,wave,int,weights,times_name,filestubname)
fid=fopen(times_name,'w+');
for s=1:numel(jd)
    fprintf(fid,'%s %.6f %.1f\n',[filestubname num2str(s) '.dat'],jd(s),weights(s));
    tempvar=cat(2,wave(:),int(s,:)');
    save([filestubname num2str(s) '.dat'],'tempvar','-ascii')
end
