function [x_nml_wave,UB] = normalize_UB( H_wave,y_wave,x,N,l,u)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
x_cplx = x(1:N)+j*x(N+1:2*N);
x_temp = x_cplx./abs(x_cplx);
for ii = 1:N
    delta = u(ii)-l(ii);
    if delta <=pi
        x_nml(ii,1) = x_temp(ii);
    else
        if ((angle(x_temp(ii))<=u(ii))&&(angle(x_temp(ii))>=l(ii)))
            x_nml(ii,1) = x_temp(ii);
        else
            d_l = abs(x_cplx(ii)-exp(j*l(ii)));
            d_u = abs(x_cplx(ii)-exp(j*u(ii)));
           if d_l <= d_u
                x_nml(ii,1) = exp(j*l(ii));
           else
                x_nml(ii,1) = exp(j*u(ii));
           end
        end
    end
end
x_nml_wave = [real(x_nml);imag(x_nml)];    
UB = objval_func(x_nml_wave,H_wave,y_wave);
end

