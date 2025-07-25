%This code sets up an expression for the dark exposure sums as a function of exposure. Modify the file names and times of exposure to get new fits.
%The file dark_exposure_fit is used in stellar_processing

files={'J8466007.fit','J8466008.fit','J8466009.fit','J8466010.fit'}%1min
%files={'J8465007.fit','J8465008.fit','J8465009.fit','J8467004.fit','J8467005.fit','J8467006.fit'}%5min
%files={'J8467008.fit','J8467009.fit'}%10min
%files={'J8467016.fit','J8467017.fit','J8467018.fit'}%12min
%files={'J8465001.fit','J8465002.fit','J8465003.fit'}%15min
%files={'J8467013.fit','J8467014.fit','J8467015.fit'}%20min
%files={'J8465004.fit','J8465005.fit','J8465006.fit','J8467010.fit','J8467011.fit','J8467012.fit'}%25min
%files={'J8464011.fit','J8464013.fit','J8464014.fit','J8464015.fit'}%30min

dark_sum=[];
for n=1:numel(files)
    eval(['matrix=raw_hercules_fitsread(''' files{n} ''');'])
    matrix1=reshape(matrix,[1,numel(matrix)]);
    matrix1=sort(matrix1);
    matrix1(end-500:end)=[];
    dark_sum(n)=sum(matrix1,'all');
end

figure
plot(dark_sum,'x')

dark_cnt(1)=mean(dark_sum(1:3));
dark_cnt(2)=mean(dark_sum(1:6));
dark_cnt(3)=mean(dark_sum(1:2));
dark_cnt(4)=mean(dark_sum(1:3));
dark_cnt(5)=mean(dark_sum(2:3));
dark_cnt(6)=mean(dark_sum(1:3));
dark_cnt(7)=mean(dark_sum(1:6));
dark_cnt(8)=mean(dark_sum(2:4));

dark_times=[1 5 10 12 15 20 25 30];
dark_times=dark_times*60;

figure
plot(dark_times,dark_cnt,'x')
[p,S,mu] = polyfit(dark_times,dark_cnt,1);
thefit=polyval(p,dark_times,S,mu);
hold on
plot(dark_times,thefit,'r')

save('dark_exposure_fit',[mu,p,S])