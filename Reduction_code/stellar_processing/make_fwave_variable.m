function [fwave1,fwave2]=make_fwave_variable(trueindbefore,trueindafter,fwave) 
%this is required for parallel processing
eval(['fwave1=fwave.th' num2str(trueindbefore(1)) ';'])
eval(['fwave2=fwave.th' num2str(trueindafter(1)) ';'])