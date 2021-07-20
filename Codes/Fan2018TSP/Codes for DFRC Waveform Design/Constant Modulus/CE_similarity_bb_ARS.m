%% Producing Fig. 8 ARS
clc;
clear all;
close all;
warning off;
N = 16;   % Antenna Number
K = 4;    % Users Number
L = 20;   % Length of Frame
power = 10^(0/10); % Total Transmit Power
amp = sqrt(power/N); % Amplitude of the Transmit Signal

% randn('state',1);
H = (randn(N,K)+j*randn(N,K))/sqrt(2); % Channel
N_pbits = 2*K;
msg_bits = randi([0,1],1,N_pbits);
% y = QPSK_mapper(msg_bits).';  % Desired Symbol
y = QPSK_mapper([0,0,0,1,1,1,1,0]).';
y_wave = sqrt(power)*[real(y);imag(y)]; % Equivalent Real Desired Symbol

for ii = 1:N
    for nn = 1:L
        X0(ii,nn) = exp(j*2*pi*ii*(nn-1)/L)*exp(j*pi*(nn-1)^2/L);  % Reference Radar Signal (LFM)
    end
end
ee = 1; % Inf Norm Similarity
H_wave = amp*[real(H),imag(H);-imag(H),real(H)]; % Equivalent Real Channel
x0 = X0(:,4); %Reference Radar Signal Vector
x0_wave = [real(x0);imag(x0)]; % Equivalent Real Radar Signal
delta = acos(1-ee^2/2);
for ii = 1:N
    l(ii,1) = angle(x0(ii))-delta;
    u(ii,1) = angle(x0(ii))+delta;     %Initialized Upper and Lower Bound
end                    
A = zeros(N,2*N);
for ii = 1:N
    A(ii,ii) = cos((l(ii)+u(ii))/2)/cos(delta);
    A(ii,ii+N) = sin((l(ii)+u(ii))/2)/cos(delta);
end                  %Initialized Linear Constraints

max_iternum = 1000; %Maximum Iteration Number
epsl = 1e-4; %Tolerence
epsl1 = 1e-6;
%-------------Parameter Initialization
[x,LB] = QCQP_LB1( H_wave,y_wave,N,l,u);          %Initialized LB and x
[x_nml1,UB1] = normalize_UB( H_wave,y_wave,x,N,l,u); %Initialized Normalization UB
[x_nml2,UB2] = QCQP_UB( H_wave,y_wave,N,l,u,x_nml1); % fmincon UB
[x_nml,UB] = QCQP_UB( H_wave,y_wave,N,l,u,x_nml2); % fmincon UB
LB_start = LB;
UB_start = UB;

prob_list = zeros(max_iternum+100,4*N+1); 
prob_list(1,:)=[x',l',u',LB];
used = 1;
lbest = LB;
ubest = UB;
x_opt = x_nml;

lb_seq = lbest;
ub_seq = ubest;

if (ubest-lbest)/abs(ubest)<epsl
    final_lb=lbest;
    final_ub=ubest;
end

iter = 2;
con = 1;
while iter<=max_iternum
    xc = prob_list(con,1:2*N)';
    lc = prob_list(con,(2*N+1):3*N )';
    uc = prob_list(con,(3*N+1):4*N)';
    x_cplx = xc(1:N)+j*xc(N+1:2*N);
    [x_nml3,~] = normalize_UB( H_wave,y_wave,xc,N,lc,uc);
    x_nml3_cplx = x_nml3(1:N)+j*x_nml3(N+1:2*N);
    x_abs = abs(x_cplx - x_nml3_cplx);
    [~,cd] = max(x_abs);
    
    
    
    
%     x_abs = abs(x_cplx);
%     [~,cd] = min(x_abs);
    xchild_left_lb=lc;
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
    tic;
    [x,lb] = QCQP_LB1( H_wave,y_wave,N,xchild_left_lb,xchild_left_ub);
    timer1(iter-1) = toc;
%     tic;
%     [xn_temp,ub_temp] = normalize_UB( H_wave,y_wave,x,N,xchild_left_lb,xchild_left_ub);
%     timer2(iter-1) = toc;
    tic;
%     [xn,ub] = QCQP_UB( H_wave,y_wave,N,xchild_left_lb,xchild_left_ub,xn_temp);
    
    [xn,ub] = normalize_UB( H_wave,y_wave,x,N,xchild_left_lb,xchild_left_ub);
    timer3(iter-1) = toc;
    
    
    
    
    if ub < ubest
       ubest=ub;
       x_opt=xn;
    end
    prob_list(used+1,:)=[x',xchild_left_lb',xchild_left_ub',lb];
    used=used+1;
    
    [x,lb] = QCQP_LB1( H_wave,y_wave,N,xchild_right_lb,xchild_right_ub);
%     [xn_temp,ub_temp] = normalize_UB( H_wave,y_wave,x,N,xchild_right_lb,xchild_right_ub);
%     [xn,ub] = QCQP_UB( H_wave,y_wave,N,xchild_right_lb,xchild_right_ub,xn_temp);
    [xn,ub] = normalize_UB( H_wave,y_wave,x,N,xchild_right_lb,xchild_right_ub);
    if ub < ubest
       ubest=ub;
       x_opt=xn;
    end
    prob_list(used+1,:)=[x',xchild_right_lb',xchild_right_ub',lb];
    used=used+1;
    
    
    [lbest,con]=min(prob_list(1:used,4*N+1));

    lb_seq(iter)=lbest;
    ub_seq(iter)=ubest;
    iter=iter+1;
    
    if ((ubest-lbest)/abs(ubest)<epsl || (ubest-lbest)<epsl1)
        final_lb=lbest;
        final_ub=ubest;
        break;
    end
    clc
    disp(['Progress - ',num2str(iter),'/',num2str(max_iternum)]); 
end
% timer_tot = sum(timer1)+sum(timer3);%+sum(timer2)
x_cplx = x_opt(1:N)+j*x_opt(N+1:2*N); %Optimal Complex Signal Vector
y_rc = amp*H.'*x_cplx; % Noise-free Received Symbol
%Constraints Check
inf_norm = norm(x_cplx-x0,Inf);
elp = abs(x_cplx);

%%
plot(1:length(lb_seq),lb_seq,'LineWidth',1.5);hold on;plot(1:length(lb_seq),ub_seq,'LineWidth',1.5);grid on;

        

            
        
        
        
            
        
