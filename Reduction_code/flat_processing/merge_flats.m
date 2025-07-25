function summed_flat = merge_flats(info,Bias)

%find all flat fields and check that they do not saturate and that they all agree with each other and then add them

%Do you want to ignore flats that are within 5000 ADU of saturation?')
ing='n';%input('(y/[n]) ->','s');

maxdatalimit=65526;

numfiles=numel(info.list);
cnt=0;
for s=1:numfiles
    eval(['temptype=info.' info.list{s} '.HERCEXPT;'])
    if strcmp(temptype,'White L')
        fprintf('reading file %s.fit\n',info.list{s})
        data = raw_hercules_fitsread([info.list{s} '.fit']);
        if strcmp(ing,'y') && max(max(data)) > maxdatalimit-5000
            continue
        else
            cnt=cnt+1;
            if exist('Bias')
                if size(Bias)==size(data)
                    data = data - Bias;
                end
            end
            if cnt==1
                summed_flat=data;
            else
                summed_flat=summed_flat+data;
            end
        end
    end
end
if ~exist('summed_flat','var')
    summed_flat=[];
end
disp('flat summing complete')