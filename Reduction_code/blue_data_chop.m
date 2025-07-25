function data = blue_data_chop(data,value)
%chop off blue orders (first 'value' pixel rows)
data(1:value,:)=[];
