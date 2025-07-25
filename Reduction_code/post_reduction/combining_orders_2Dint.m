function [fullwave,fullint_rel] = combining_orders_2Dint(data,step,ordmin,ordmax)
%make final wavelength axis
eval(['minwave1=min(data.wave' num2str(ordmin) ');'])
eval(['minwave2=min(data.wave' num2str(ordmax) ');'])

if minwave1>minwave2
    %then redo order numbering so increasing order number corresponds to
    %increasing wavelength
    cnting=0;
    for J=ordmin:ordmax
        cnting=cnting+1;
        NEW=ordmax+1-cnting;
        eval(['TMP.wave' num2str(NEW) '=data.wave' num2str(J) ';'])
        eval(['TMP.int' num2str(NEW) '=data.int' num2str(J) ';'])
        eval(['TMP.weights' num2str(NEW) '=data.weights' num2str(J) ';'])
    end
    data=TMP;
    clear TMP
end

eval(['fullwave=rounding( min(data.wave' num2str(ordmin) '), step) : step : rounding( max(data.wave' num2str(ordmax) '), step);'])
eval(['numobs=size(data.int' num2str(ordmin) ',1);'])
eval(['fullint_rel=zeros( numel(data.int' num2str(ordmin) '(:,1)), length(fullwave) );'])

for OBS=1:numobs
    fullweight=zeros( 1, length(fullwave) );
    weights_even=fullweight;
    weights_odd=fullweight;
    int_even=fullweight;
    int_odd=fullweight;

    for s=ordmin:ordmax
        eval(['wave=data.wave' num2str(s) ';'])
        eval(['int=data.int' num2str(s) '(OBS,:);'])
        eval(['weight=data.weights' num2str(s) '(OBS,:);'])
        
        S = rounding(fullwave,step) >= min(rounding(wave,step)) & rounding(fullwave,step) <= max(rounding(wave,step));
        if sum(S)~=numel(wave)
            S = rounding(fullwave,step) >= ceil(min(wave)/step)*step & rounding(fullwave,step) <= ceil(max(wave)/step)*step;
        end
        fullweight(S)=fullweight(S)+weight;
        if round(s)/2==round(s/2)
            weights_even(S)=weight;
            int_even(S)=int;
        else
            weights_odd(S)=weight;
            int_odd(S)=int;
        end
    end
    fullint_rel(OBS,:)=int_even.*weights_even./fullweight + int_odd.*weights_odd./fullweight;
    clear fullweight weights_even weights_odd int_even int_odd
end

function output = rounding(input,tonearest)
output = round(input/tonearest) * tonearest;

function nearest = findnearestlowest(vector,value)
test=vector-value;
nearest=find( abs(test) == min(abs(test)) );
if numel(nearest) > 1
    nearest=min(nearest);
end

function nearest = findnearesthighest(vector,value)
test=vector-value;
nearest=find( abs(test) == min(abs(test)) );
if numel(nearest) > 1
    nearest=max(nearest);
end
