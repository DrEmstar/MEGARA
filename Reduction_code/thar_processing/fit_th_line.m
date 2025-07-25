function [x,resnorm,residual,exitflag,fitted] = fit_th_line(xdata,ydata,lb,ub,startpoint)

%fit simple Gaussian to data
fun = @(x,xdata)x(1)*exp(-(((xdata-x(2))*(2*sqrt(log(2)))/x(3)).^2))+x(4)+x(5)*xdata;%,'x','xdata';
options = optimset('Display','off');
[x,resnorm,residual,exitflag] = lsqcurvefit(fun,startpoint,xdata,ydata,lb,ub,options);
fitted=fun(x,xdata);