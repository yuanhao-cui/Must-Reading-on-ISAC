function [ PD ] = PD_Orthogonal( X,a,SNR )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[N,L] = size(X);
Rs = X*X'/L;
rou = SNR*abs(a'*Rs.'*a)^2;
Pfa = 1e-7;
ita = chi2inv(1-Pfa,2);
PD = 1-ncx2cdf(ita,2,rou);
end

