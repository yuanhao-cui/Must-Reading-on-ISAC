function [ x, obj_val,exitflag ] = QCQP_LB2( H_wave,y_wave,N,l,u,obj_prev)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
Q = 2*H_wave*H_wave.';
f = -2*H_wave*y_wave;
c = y_wave.'*y_wave;
A = zeros(N,2*N);
for ii = 1:N
    A(ii,ii) = cos((l(ii)+u(ii))/2);
    A(ii,ii+N) = sin((l(ii)+u(ii))/2);         %Linear Constraint
end                         


function [y,grady] = quadobj(x,Q,f,c)         %Objective Function
    y = 1/2*x.'*Q*x+f.'*x+c;
    if nargout>1
        grady = Q*x + f;
    end
end

function [y_con,yeq,grady,gradyeq] = quadconstr(x)      %Nonlinear Constraint
    for nn = 1:N
        y_con(nn) = x(nn)^2+x(nn+N)^2 - 1;
    end
    y_con(N+1) = 1/2*x.'*Q*x+f.'*x+c-obj_prev;
    yeq = [];
    if nargout>2
       gradyeq = zeros(length(x),N+1);
       for nn = 1:N
           grady(nn,nn) = 2*x(nn);
           grady(nn+N,nn) = 2*x(nn+N);
       end
       grady(:,N+1) = Q*x + f;
    end
    gradyeq = [];
end

Acon = -A;

for ii = 1:N
    b(ii) = -cos((u(ii)-l(ii))/2);
end


% options = optimoptions(@fmincon,'Algorithm','sqp',...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
%     'Display','off','OptimalityTolerance',1e-6,'MaxIterations',400);

options = optimoptions('fmincon','Algorithm','sqp', 'GradObj','on');

fun = @(x)quadobj(x,Q,f,c);
nonlconstr = @(x)quadconstr(x);
x1 = ones(2*N,1);
[x,ignore,exitflag] = fmincon(fun,x1,Acon,b,[],[],[],[],nonlconstr,options);
obj_val = objval_func(x,H_wave,y_wave);
end

