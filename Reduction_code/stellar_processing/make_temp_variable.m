function [temp,jdexp] = make_temp_variable(info,s)
%this is required by the parallel loop
eval(['temp=info.' info.list{s} '.HERCEXPT;'])
eval(['jdexp=info.' info.list{s} '.MJD_OBS;'])