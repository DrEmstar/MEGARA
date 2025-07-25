function [line] = get_the_line_profile(star,starname,type,s,files)
A = fitsread([files.list{s} '.fit']); %load data
%check_the_file(A,s) %check quality of file
y=[];
try
    %creat matrix with 80 rows to average around 2000 (middle)  
    y = [y;A(2000,:)]; 
catch
    error(strcat('There is a size issue with file: ',files.list{s} ,'.fit'))
end 
for n=1960:2039
    %creat matrix with 80 rows to average around 2000 (middle)
    y = [y;A(n,:)]; 
end

avrg = [];
for n=1:4130
    %AVERAGE ALONG COLOUMNS
    p =mean(y(:,n));
    avrg=cat(2,avrg,p);%average through the coloums
end
%invert array and make for y
y=avrg';
%index the double array for use later           
y = [double(y(:,1))]; 
%get data to lie on the floor
y=y-min(y);
%normalise the data
line = y/max(y);
end