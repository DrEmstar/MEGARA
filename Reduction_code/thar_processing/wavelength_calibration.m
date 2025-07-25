function wavfit=wavecal(thfitinfo,order,air)

c = 299792458; %speed of light (km/s)
chosen.orders=order(thfitinfo(:,1));
chosen.mlam=order(thfitinfo(:,1)).*air(thfitinfo(:,1));
chosen.xpos=thfitinfo(:,4);
orig.orders=chosen.orders;
orig.mlam=chosen.mlam;
orig.xpos=chosen.xpos;

x=[];y=[];z=[];
xcen=min(chosen.orders);ycen=min(chosen.xpos);zcen=min(chosen.mlam);
xscale=max(abs(chosen.orders));yscale=max(abs(chosen.xpos));zscale=max(abs(chosen.mlam));
x=(chosen.orders-xcen)/xscale;y=(chosen.xpos-ycen)/yscale;z=(chosen.mlam-zcen)/zscale;

while 1
    p = polynomial_fit_2d(x,y,z,6,4);
    zfit = polynomial_val_2d(p,x,y);
    ZFIT=zfit*zscale+zcen;
    resid=chosen.mlam-ZFIT;
    S=abs(resid)==max(abs(resid));
    er=std(resid);
    rver=er/mean(chosen.mlam)*c;
    if  rver<= 100 %m/s
        fprintf('%.0f lines chosen, %.0f lines rejected\n\n',numel(ZFIT),numel(orig.mlam)-numel(ZFIT))
        wavfit.p=p;
        wavfit.xcen=xcen;
        wavfit.ycen=ycen;
        wavfit.zcen=zcen;
        wavfit.xscale=xscale;
        wavfit.yscale=yscale;
        wavfit.zscale=zscale;
        break
    end
    chosen.xpos(S)=[];chosen.orders(S)=[];chosen.mlam(S)=[];
    x(S)=[];y(S)=[];z(S)=[];
end