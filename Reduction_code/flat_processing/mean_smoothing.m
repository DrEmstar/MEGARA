function out=mean_smoothing(x,boxsize)
out=zeros(size(x));
for s=1:numel(x)
    rang=s-boxsize:s+boxsize;
    rang(rang<1)=[];
    rang(rang>numel(x))=[];
    out(s)=mean(x(rang));
end
