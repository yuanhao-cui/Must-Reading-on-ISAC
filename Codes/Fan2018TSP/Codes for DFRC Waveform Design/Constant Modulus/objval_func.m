function [ objval ] = objval_func(x,H_wave,y_wave)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
objval = norm(H_wave.'*x-y_wave,2)^2;
end

