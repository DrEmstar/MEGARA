function res = testing_deconvolution(a,b)
 s1=fft(a-mean(a),2^14);
 s2=fft(b-mean(b),2^14);

X=s1./(abs(s2/max(abs(s2))));

if X(1)==Inf || X(1)==-Inf
    X(1)=0;
end

XX1=(real(s1).*real(s2)+imag(s1).*imag(s2))./(real(s2).^2+imag(s2).^2);
XX2=(imag(s1).*real(s2)-real(s1).*imag(s2))./(real(s2).^2+imag(s2).^2);
XX=complex(XX1,XX2);

res=real(ifft(X));
res=res(1:numel(a));
res=res+mean(a);