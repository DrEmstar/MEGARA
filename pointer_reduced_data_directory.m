%to be run in the Jfiles_bydate folder

d = dir;
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
numfold=numel(nameFolds);
for n=1:numfold
    foldnam=nameFolds{n};
    eval(['cd ' foldnam(1:3) ''' ' foldnam(5:end) ''''])
    try
        dd=dir('J*.mat');
        for m=1:numel(dd)
            movefile(dd(m).name,'../../star_files/')
        end
        
    catch
        
    end
    cd ../
end
cd ../star_files/
try
    ddd=dir('J*prc.mat');
    for l=1:numel(ddd)
        movefile(ddd(l).name,'processing_files')
    end
    
catch
end

