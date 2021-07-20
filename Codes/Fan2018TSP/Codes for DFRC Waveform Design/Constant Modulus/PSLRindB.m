function [ PSLR ] = PSLRindB( x0 )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
f = pc_freqwin( x0 );
pks = findpeaks(f);
pks = sort(pks);
PSLR = pks(length(pks))-pks(length(pks)-1);
end

