function [centre, widths, sse, F] = fitlingauss(xdata, ydata, start_point)

%this function is based on the "Fitting a Curve to Data" help topic in the
%MATLAB help

% Call fminsearch with a fixed starting point.
% this reduction version always has a specified width startpoint ...
if ~exist('start_point','var')
    start_point = [max(ydata), mean(xdata), 0, 0, 5]; %assuming centred and scaled then these should be reasonable starting points
end

%parameters are:
%[amp xshift slope yshift width];

model = @lingaussfun;
options=optimset('Display','off');
estimates = fminsearch(model, start_point, options);

sse = lingaussfun(estimates);

F = estimates(1)*exp( -( ((xdata-estimates(2))*2*sqrt(log(2))/estimates(5)).^2 ) )+estimates(3)*xdata+estimates(4);

centre=estimates(2);
widths=estimates(5);

    function sse = lingaussfun(params)%[sse, FittedCurve] = lingaussfun(params)
        a = params(1);b = params(2);c = params(3);d = params(4);e = params(5);
        FittedCurve = a*exp( -( ((xdata-b)*2*sqrt(log(2))/e).^2 ) )+c*xdata+d; %the best fit curve shape - in this case a gaussian on a sloped line
        ErrorVector = FittedCurve - ydata; %least-squares error
        sse = sum(ErrorVector .^ 2); %used for the fminsearch optimisation
    end
end
