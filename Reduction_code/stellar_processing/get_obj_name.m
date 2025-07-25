function obj_name = get_obj_name(filename)
header = get_header(filename);
for s=1:numel(header)
    if ~isempty(strfind(header{s}, 'OBJECT'))
        try
            obj_name=header{s};
            try
                obj_name=obj_name(10:31); %trim off extras
            catch
                obj_name=obj_name(10:end); %trim off extras
            end

            obj_name(double(obj_name)==32)=[]; %trim off spaces
            obj_name(double(obj_name)==39)=[]; %trim off apostrophes
        catch
            disp('no OBJECT header found!')
        end
    end
end
