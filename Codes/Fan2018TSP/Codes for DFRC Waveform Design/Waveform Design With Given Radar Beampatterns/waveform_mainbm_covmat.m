function [ R ] = waveform_mainbm_covmat( N,a,theta_target,theta,power )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
delta=theta(2)-theta(1);
l=ceil((theta_target+pi/2*ones(1,5))/delta+ones(1,5));
cvx_solver sedumi
cvx_begin quiet
variable t
variable R(N,N) hermitian semidefinite
minimize -t
subject to
for i=1:l(1)-1
    real(a(:,l(3))'*R*a(:,l(3))-a(:,i)'*R*a(:,i))>=t;
end
for ii=l(5)+1:length(theta)
    real(a(:,l(3))'*R*a(:,l(3))-a(:,ii)'*R*a(:,ii))>=t;
end
for i=l(3):l(5)-1
    real(a(:,i)'*R*a(:,i)-a(:,i+1)'*R*a(:,i+1))>=0.001*i;
end
for i=l(1):l(3)-1
    real(a(:,i+1)'*R*a(:,i+1)-a(:,i)'*R*a(:,i))>=0.001*i;
end
a(:,l(2))'*R*a(:,l(2)) == 0.5*a(:,l(3))'*R*a(:,l(3));
a(:,l(4))'*R*a(:,l(4)) == 0.5*a(:,l(3))'*R*a(:,l(3));
diag(R) == ones(N,1)*power/N;
% trace(R)==power;
cvx_end
end

