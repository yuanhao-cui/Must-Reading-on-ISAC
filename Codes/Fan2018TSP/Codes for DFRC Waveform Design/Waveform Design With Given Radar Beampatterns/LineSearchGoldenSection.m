function [t] = LineSearchGoldenSection(func,LB,UB,EPSILON)
%
% Golden Section search for an approximate minimum of unimodal 
%         function f(t) using only function evaluations
%
% length of location interval is reduced by 0.318 each iteration 
%

narginchk(3,4);
if nargin == 3
    EPSILON = max(abs(UB) + abs(LB),2)*0.5e-4;
end
GR = (sqrt(5)-1)/2; % 0.618

l = UB - (UB-LB)*GR;
r = LB + (UB-LB)*GR;
fl = func(l);
fr = func(r);
while UB - LB >= EPSILON
    if fl < fr
        UB = r;
        r = l;
        fr = fl;
        l = UB - (UB-LB)*GR;
        fl = func(l);
    else
        LB = l;
        l = r;
        fl = fr;
        r = LB + (UB-LB)*GR;
        fr = func(r);
    end
end

t = (UB+LB)/2;

end