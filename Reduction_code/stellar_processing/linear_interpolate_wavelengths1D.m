function wavelength = linear_interpolate_wavelengths1D(fwave1,fwave2,jd1,jd2,jdexp,m)

for s=1:numel(m)
    eval(['wave1=fwave1.order_' num2str(m(s)) ';'])
    eval(['wave2=fwave2.order_' num2str(m(s)) ';'])
    wave1=wave1(:);
    wave2=wave2(:);
    if numel(wave1)~=numel(wave2)
        disp('different number of elements in each order - you should use the same flat-field for reduction to avoid this!')
        return
    end
    if jd2-jd1 ~= 0
        P = (jdexp-jd1)/(jd2-jd1);
    else
        P = 0;
    end
    eval(['wavelength.order_' num2str(m(s)) '=wave1 + P*(wave2-wave1);'])
end

