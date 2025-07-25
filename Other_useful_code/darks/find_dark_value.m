function dark=find_dark_value(exp_time)

load('dark_exposure_fit.mat','p','S','mu')
dark=polyval(p,exp_time,S,mu);