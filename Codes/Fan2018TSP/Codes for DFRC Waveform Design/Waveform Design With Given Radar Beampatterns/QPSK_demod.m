function [ symbol_demod ] = QPSK_demod_mat(y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[M,N]=size(y);
for ii = 1:M
    for jj = 1:N
        if  real(y(ii,jj))>0
            if imag(y(ii,jj))>0
                symbol_demod(ii,jj)=(1+j)/sqrt(2);
            else
                symbol_demod(ii,jj)=(1-j)/sqrt(2);
            end
        else
            if imag(y(ii,jj))>0
                symbol_demod(ii,jj)=(-1+j)/sqrt(2);
            else
                symbol_demod(ii,jj)=(-1-j)/sqrt(2);
            end
        end
    end
end

