function [ R ] = waveform_design_multibm_covmat( Pd_theta,N,L,a,theta,power)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cvx_solver sedumi
cvx_begin quiet
variable bt
variable R(N,N) hermitian semidefinite
expression u(length(theta))
for i=1:length(theta)
    u(i)=(bt*Pd_theta(i)-a(:,i)'*R*a(:,i));
end
minimize norm(u,2)
subject to
% diag(R)==ones(N,1)*power/N;
trace(R)==power;
cvx_end
end

