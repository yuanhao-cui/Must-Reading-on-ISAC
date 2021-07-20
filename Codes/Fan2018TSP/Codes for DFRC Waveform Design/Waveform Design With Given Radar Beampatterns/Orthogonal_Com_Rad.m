function [ X ] = Orthogonal_Com_Rad( H,Y,power )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[N,K] = size(H);
[~,L] = size(Y);
u = sqrt(power/N);
[U,S,V] = svd(sqrt(L/N)*power*conj(H)*Y);
X = u*sqrt(L)*U*eye(N,L)*V';
end

