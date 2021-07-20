function [ X0 ] = chirp_generator( N,fc,Tp,B )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
u = B/Tp;
fs = 4*B;
sp_num = fs*Tp;
t = linspace(0,Tp,sp_num);
for ii = 1:N
     X0(ii,:) = exp(j*2*pi*10*ii*fc*(t-1)).*exp(j*pi*u*(t-1).^2);
end
end

