function JD = convert_date_to_JD(year,month,day,time)
%only for obs taken after the start of the year 2000!
%time in ??:??:??.?? format entered as a string - MUST be UT 24 hour time
%year month day entered as numbers

%The Julian date for CE  2000 January  1st  00:00:00.0 UT is
%JD 2451544.50000

%days to the start of the year of obs
if year < 100 %two digit year was entered
    year = year+2000;
end
yearspast = year-2000;
yeardays=365*yearspast+ floor((yearspast-1)/4) + 1;

allmonthdays = [31 28 31 30 31 30 31 31 30 31 30 31];
if month > 1
    simplemonthdays = sum(allmonthdays(1:month-1));
else
    simplemonthdays = 0;
end
L = isleapyear(year);
if month > 2 && L == 1
    monthdays = simplemonthdays + 1;
else
    monthdays = simplemonthdays;
end
totaldays=yeardays+monthdays+day-1; %the day today hasn't finished yet so -1

hours = str2num(time(1:2));
minutes = str2num(time(4:5));
seconds = str2num(time(7:11));
if hours > 0
    hoursecs = hours*3600;
else
    hoursecs = 0;
end
if minutes > 0
    minutesecs = minutes*60;
else
    minutesecs = 0;
end
totalsecs = hoursecs+minutesecs+seconds;
dayfraction = totalsecs/86400; % 86400 seconds in a day close approx.

JD = 2451544.5 + totaldays + dayfraction;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = isleapyear(year)
%only for obs taken after the start of the year 2000!
if round(year/4) == year/4
    L = 1;
else
    L = 0;
end
