function [x_cplx,opt_val] = CE_similarity_ComRad_benchmark( H,y,power,ee,x0,cle )
% Branch and Bound Method for ComRad CEP
%   Detailed explanation goes here
N = length(x0);
amp = sqrt(power/N); % Amplitude of the Transmit Signal
y_wave = sqrt(power*cle)*[real(y);imag(y)]; % Equivalent Real Desired Symbol
H_wave = amp*[real(H),imag(H);-imag(H),real(H)]; % Equivalent Real Channel
x0_wave = [real(x0);imag(x0)];
delta = acos(1-ee^2/2);
for ii = 1:N
    l(ii,1) = angle(x0(ii))-delta;
    u(ii,1) = angle(x0(ii))+delta;     %Initialized Upper and Lower Bound
end           
l0 = l;
u0 = u;

A = zeros(N,2*N);
for ii = 1:N
    A(ii,ii) = cos((l(ii)+u(ii))/2)/cos(delta);
    A(ii,ii+N) = sin((l(ii)+u(ii))/2)/cos(delta);
end     

max_iternum = 200; %Maximum Iteration Number
epsl = 1e-3; %Tolerence
epsl1 = 1e-6;
%-------------Parameter Initialization
[x,obj_val] = QCQP_LB1( H_wave,y_wave,N,l,u);          %Initialized LB and x
iter = 1;

while iter<=max_iternum
    tr = (l+u)/2;
    x_cplx = x(1:N)+j*x(N+1:2*N);
    obj_prev = obj_val;
    obj(iter) = obj_val;
    for ii = 1:N
        theta = angle(x_cplx(ii));
        if theta>=tr(ii)
            l(ii) = tr(ii);
        else
            u(ii) = tr(ii);
        end
    end
    [x,obj_val,exitflag] = QCQP_LB2( H_wave,y_wave,N,l,u,obj_prev);
    iter = iter+1;
    if exitflag == -2
        x_opt = [real(x_cplx);imag(x_cplx)];
        break;
    elseif (abs(obj_prev/obj_val-1)<=epsl)||(abs(obj_prev-obj_val)<=epsl1)
        x_opt = x;
        break;
    end
end

x_opt1 = normalize_UB( H_wave,y_wave,x_opt,N,l,u);
x_cplx = x_opt1(1:N)+j*x_opt1(N+1:2*N);
opt_val = objval_func(x_opt1,H_wave,y_wave);



end

