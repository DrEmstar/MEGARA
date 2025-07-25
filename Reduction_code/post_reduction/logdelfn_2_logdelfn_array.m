function ldelfn_array = logdelfn_2_logdelfn_array(ldelfn,logwavelength)

ldelfn_array=zeros(size(logwavelength));
for s=1:numel(ldelfn(:,1))
    T=find_nearest(logwavelength,ldelfn(s,1));
    ldelfn_array(T)=ldelfn_array(T)+ldelfn(s,2);
end

function nearest = find_nearest(vector,value)
test=vector-value;
nearest=find( abs(test) == min(abs(test)) );
if numel(nearest) > 1
    nearest=nearest(1); 
end

