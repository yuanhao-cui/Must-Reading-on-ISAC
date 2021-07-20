function [x_cplx,opt_val,lb_seq,ub_seq] = CE_similarity_ComRad( H,y,power,ee,x0,cle )
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

A = zeros(N,2*N);
for ii = 1:N
    A(ii,ii) = cos((l(ii)+u(ii))/2)/cos(delta);
    A(ii,ii+N) = sin((l(ii)+u(ii))/2)/cos(delta);
end     

max_iternum = 200; %Maximum Iteration Number
epsl = 1e-3; %Tolerence
epsl1 = 1e-6;

[x,LB] = QCQP_LB1( H_wave,y_wave,N,l,u);          %Initialized LB and x
[x_nml1,UB1] = normalize_UB( H_wave,y_wave,x,N,l,u); %Initialized Normalization UB
[x_nml2,UB2] = QCQP_UB( H_wave,y_wave,N,l,u,x_nml1); % fmincon UB
[x_nml,UB] = QCQP_UB( H_wave,y_wave,N,l,u,x_nml2); % fmincon UB
LB_start = LB;
UB_start = UB;

prob_list = zeros(max_iternum+100,4*N+1);          %Problem list initialization
prob_list(1,:)=[x',l',u',LB];
used = 1;
lbest = LB;
ubest = UB;
x_opt = x_nml;

lb_seq = lbest;
ub_seq = ubest;

% if (ubest-lbest)/abs(ubest)<epsl
%     final_lb=lbest;
%     final_ub=ubest;
% end

iter = 2;
con = 1;


while iter<=max_iternum
    xc = prob_list(con,1:2*N)';                                   % Pick a problem having the smallest LB
    lc = prob_list(con,(2*N+1):3*N )';
    uc = prob_list(con,(3*N+1):4*N)';
    x_cplx = x(1:N)+j*x(N+1:2*N);               
    %     x_abs = abs(x_cplx);
%     [x_abs_min,cd] = min(x_abs);
    
    [x_nml3,~] = normalize_UB( H_wave,y_wave,x,N,lc,uc);          % Calculate the UB of the problem
    x_nml3_cplx = x_nml3(1:N)+j*x_nml3(N+1:2*N);
    x_abs = abs(x_cplx - x_nml3_cplx);                            % Branching from the x(n) having the largest gap between x_u and x_l
    [~,cd] = max(x_abs);
    
    xchild_left_lb=lc;                                            % Generate two sub-problems (left child and right child) at the chosen x(n)
    xchild_left_ub=uc;
    xchild_right_lb=lc;
    xchild_right_ub=uc;
    tr=(lc(cd)+uc(cd))/2;
    xchild_left_ub(cd)=tr;
    xchild_right_lb(cd)=tr;
    
    if con < used
        prob_list(con,:) = prob_list(used,:);
        used=used-1;
    else
        used=used-1;
    end
    [x,lb] = QCQP_LB1( H_wave,y_wave,N,xchild_left_lb,xchild_left_ub);        % Compute the LB and the associated solution of the left sub-problem

    
%     [xn_temp,ub_temp] = normalize_UB( H_wave,y_wave,x,N,xchild_left_lb,xchild_left_ub);
%     [xn,ub] = QCQP_UB( H_wave,y_wave,N,xchild_left_lb,xchild_left_ub,xn_temp);
    [xn,ub] = normalize_UB( H_wave,y_wave,x,N,xchild_left_lb,xchild_left_ub); % Compute the UB and the associated solution of the left sub-problem
    
    
    
    
    if ub < ubest                                                   % If the UB is lower than the current ubest, replace ubest with UB
       ubest=ub;
       x_opt=xn;
    end
    prob_list(used+1,:)=[x',xchild_left_lb',xchild_left_ub',lb];              % Insert the associated LB and solutions into the problem list
    used=used+1;
    
    [x,lb] = QCQP_LB1( H_wave,y_wave,N,xchild_right_lb,xchild_right_ub);                         % Compute the LB and the associated solution of the right sub-problem
%     [xn_temp,ub_temp] = normalize_UB( H_wave,y_wave,x,N,xchild_right_lb,xchild_right_ub);
%     [xn,ub] = QCQP_UB( H_wave,y_wave,N,xchild_right_lb,xchild_right_ub,xn_temp);
    [xn,ub] = normalize_UB( H_wave,y_wave,x,N,xchild_right_lb,xchild_right_ub); % Compute the UB and the associated solution of the right sub-problem
    if ub < ubest                    % If the UB is lower than the current ubest, replace ubest with UB
       ubest=ub;
       x_opt=xn;
    end
    prob_list(used+1,:)=[x',xchild_right_lb',xchild_right_ub',lb]; % Insert the associated LB and solutions into the problem list
    used=used+1;
    
    
    [lbest,con]=min(prob_list(1:used,4*N+1)); % Replace lbest with the smallest LB in the list, mark its index in the problem list as con

    lb_seq(iter)=lbest;           %Generate LB and UB sequences
    ub_seq(iter)=ubest;
    iter=iter+1;
    
    if ((ubest-lbest)/abs(ubest)<epsl || (ubest-lbest)<epsl1)               %Convergence condition
        final_lb=lbest;
        final_ub=ubest;
        break;
    end
%     clc
%     disp(['Progress - ',num2str(iter),'/',num2str(max_iternum)]); 
end
x_cplx = x_opt(1:N)+j*x_opt(N+1:2*N);          %Optimal x

opt_val = objval_func(x_opt,H_wave,y_wave);  %Optimal objective function



end

