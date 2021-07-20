function [ f ] = pc_freqwin( x0 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
w = taylorwin(length(x0)).';
s = x0.*w;
f1 = abs(fftshift((ifft(fft(s).*conj(fft(s))))));
f = 20*log10(f1/max(f1));
end

