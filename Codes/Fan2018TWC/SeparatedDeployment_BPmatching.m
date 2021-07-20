% Simulating the Separated Deployment beamformers presented in 
%   "MU-MIMO Communications with MIMO Radar: From Co-existence to Joint
%   Transmission" - Liu et. al. 2018

% Nate Raymondi, 8/13/2020

% Set up the system -------------------------------------------------------
Nc = 16;             % number of comms Tx antennas
Nr = 16;             % number of radar Tx antennas
fc = 5e9;           % assumed center freq
c = 3e8;
lamb = c/fc;
spacing = lamb/2;
radAntLoc = spacing*[0:Nr-1];
commAntLoc = spacing*[0:Nc-1];
Pc = 5;             % comms power budget
Pr = 5;             % rad power budget
N0 = 1;             % noise variance


% create ALL steering vectors ---------------------------------------------
M = 181;                                    % number of points on the angular grid
angleSpace = linspace(-pi/2, pi/2, M);
angleSpaceDeg = linspace(-90, 90, M);

% rad ---------------------------------------------------------------------
a1 = zeros(Nr, length(angleSpace));         % radar antenna steering vecs
for j = 1:Nr
    a1(j,:) = exp((1i * 2 * pi * radAntLoc(j) / lamb) .* sin(angleSpace));
end

% get the actual target steering vecs for convenience
K = 2;              % number of radar targets
radAngles = [25 50];          
aRad = zeros(Nr,K); R = zeros(Nr,Nr,K); alpha = zeros(1,K);
for k = 1:K
    aRad(:,k) = a1(:,90+radAngles(k));
    R(:,:,k) = aRad(:,k) * aRad(:,k)';
    %alpha(k) = rand + 1i*rand;             % target and channel coeffs if necessary
end

% comm --------------------------------------------------------------------
a2 = zeros(Nc, length(angleSpace));         % comm antenna steering vecs
for j = 1:Nc
    a2(j,:) = exp((1i * 2 * pi * commAntLoc(j) / lamb) .* sin(angleSpace));
end

% get the actual user steering vecs for convenience 
N = 2;              % number of comms users
commAngles = [90 -15];   
Gamma = Pc*[0.5 0.2];
g = zeros(Nc,N); G = zeros(Nc,Nc,N);            % channel from comms antennas to users
f = zeros(Nr,N); F = zeros(Nr,Nr,N);            % channel from rad antennas to users
for n = 1:N
    g(:,n) = a2(:,90+commAngles(n));
    G(:,:,n) = g(:,n) * g(:,n)';
    
    f(:,n) = a1(:,90+commAngles(n));
    F(:,:,n) = f(:,n) * f(:,n)';
end

% introduce some further notation -----------------------------------------
a = [a1; a2];                   % size Nr+Nc x M
% NOTE: the way we have contructed our steering vectors and stacked them
% into the matrices a, a1, a2 - these are equal to the A, A1, A2 from Eq.17



% FIRST METHOD - BEAMPATTERN DESIGN ---------------------------------------
% create our desired beampattern to approximate 
Pdesired = zeros(length(angleSpace),1); index = 1;
beamwidth = 10;
for i = -90:90
    if min(abs(i-radAngles(:))) <= beamwidth/2
        Pdesired(index) = 10;
    else
        Pdesired(index) = 0;
    end
    index = index+1;
end
%figure; plot(angleSpaceDeg,Pdesired); title('Desired Beampattern')


% solve the optimization in Eq.12
cvx_begin quiet
    variable R1(Nr,Nr) hermitian
    variable alpha 
    minimize( sum_square_abs(alpha.*Pdesired - diag(a1' * R1 * a1)) );
    subject to
        diag(R1) == Pr*ones(Nr,1)/Nr;
        alpha >= 0;
        for n = 1:N
            trace( conj(f(:,n)) * transpose(f(:,n)) * R1 ) == 0;
        end
    R1 == semidefinite(Nr);
cvx_end
% NOTE: I am fairly certain the expression of Eq.12 in the paper is
% incorrect. The 'a' in the objective function should be 'a1'


% substitute R1 into Eq.19 and solve
%   must be solved jointly since Wi depends on all other Wj
cvx_begin quiet
    variable W1opt(Nc,Nc) hermitian
    variable W2opt(Nc,Nc) hermitian
    variable sigma
    minimize( square_pos(norm( diag( a2'*(W1opt+W2opt)*a2 - sigma*(a1'*R1*a1) ) ,2)) );
    subject to
        % SINR constraints for each user
        %abs(trace( conj(g(:,1))*transpose(g(:,1))*W1 ))/abs(( trace( conj(g(:,2))*transpose(g(:,2))*W1 ) + ...
           %trace( conj(f(:,1))*transpose(f(:,1))*R1 ) + N0 )) >= Gamma(1);
        real(trace( conj(g(:,1))*transpose(g(:,1))*W1opt ) - Gamma(1)*( trace( conj(g(:,2))*transpose(g(:,2))*W1opt ) + ...
           trace( conj(f(:,1))*transpose(f(:,1))*R1 ) )) >= N0*Gamma(1);
        real(trace( conj(g(:,2))*transpose(g(:,2))*W2opt ) - Gamma(2)*( trace( conj(g(:,1))*transpose(g(:,1))*W2opt ) + ...
           trace( conj(f(:,2))*transpose(f(:,2))*R1 ) )) >= N0*Gamma(2);
        % Power constraint
        trace(W1opt + W2opt) <= Pc;
        sigma >= 0;
    W1opt == semidefinite(Nc);
    W2opt == semidefinite(Nc);
cvx_end
% NOTE - we are solving the SDR for W1 and W2, so we will have to recover
% w1 and w2


% Do randomization to get our feasible beamforming vectors w1 and w2-------
nRand = 10;
w1 = zeros(Nc,nRand); w2 = zeros(Nc,nRand);
W1cov = zeros(Nc,Nc,nRand); W2cov = zeros(Nc,Nc,nRand); 
clear u1F; clear u2F; 

% generate and scale to meet power constraints
for L = 1:nRand
    % generate nRand
    w1(:,L) = mvnrnd(zeros(Nc,1),W1opt) + 1i*mvnrnd(zeros(Nc,1),W1opt);
    w2(:,L) = mvnrnd(zeros(Nc,1),W2opt) + 1i*mvnrnd(zeros(Nc,1),W2opt);

    % scale them so they all adhere to norm constraint
    w1(:,L) = sqrt(trace(W1opt))*w1(:,L)/sqrt(w1(:,L)'*w1(:,L));
    w2(:,L) = sqrt(trace(W2opt))*w2(:,L)/sqrt(w2(:,L)'*w2(:,L));

    W1cov(:,:,L) = w1(:,L)*w1(:,L)';
    W2cov(:,:,L) = w2(:,L)*w2(:,L)';
end

% check all pairs to see if user 1 SINR constraints are met
index1 = 1; index2 = 1;
for i = 1:nRand
    for j = 1:nRand
        if real(trace( conj(g(:,1))*transpose(g(:,1))*W1cov(:,:,i) ) - Gamma(1)*( trace( conj(g(:,2))*transpose(g(:,2))*W1cov(:,:,i) ) + ...
           trace( conj(f(:,1))*transpose(f(:,1))*R1 ) )) >= N0*Gamma(1)
            % save this pair as feasible for user 1
            u1F(:,1,index1) = w1(:,i); u1F(:,2,index1) = w2(:,j);
            index1 = index1 + 1;
        end
    end
end

% check all user 1 feasible pairs to see if they're user 2 feasible 
for i = 1:size(u1F,3)
    if real(trace( conj(g(:,2))*transpose(g(:,2))*(u1F(:,2,i)*u1F(:,2,i)') ) - Gamma(2)*( trace( conj(g(:,2))*transpose(g(:,2))*(u1F(:,1,i)*u1F(:,1,i)') ) + ...
           trace( conj(f(:,2))*transpose(f(:,2))*R1 ) )) >= N0*Gamma(2)
       % save this pair as feasible for user 2
        u2F(:,1,index2) = u1F(:,1,i); u2F(:,2,index2) = u1F(:,2,i);
        index2 = index2 + 1;
    end
end
        
% find how many feasible sets we have, choose the one with the smallest
% objective function value
numFeasible = size(u2F,3);           % this number better be [0,nRand^N] lol. Yes there are repeats but who cares
WcovSum = zeros(Nc,Nc,numFeasible); obj = zeros(1,numFeasible);
for i = 1:size(u2F,3)
    WcovSum(:,:,i) = (u2F(:,1,i)*u2F(:,1,i)') + (u2F(:,2,i)*u2F(:,2,i)');
    obj(i) = ( norm( diag( a2'*WcovSum(:,:,i)*a2 - sigma*a1'*R1*a1 ) ) )^2;
end
[minVal,minIndex] = min(obj);
w1 = u2F(:,1,minIndex); W1 = (w1*w1');
w2 = u2F(:,2,minIndex); W2 = (w2*w2');
Wsum = W1 + W2;


% the overall transmitted covariance Eq.14 --------------------------------
C = blkdiag(R1,Wsum);           % should be size (Nr+Nc) x (Nr+Nc)

% overall transmitted beampattern gain Eq.15
Pd = zeros(size(angleSpace)); BPrad = zeros(size(angleSpace)); BPcomm = zeros(size(angleSpace));
for i = 1:length(angleSpace)
    Pd(i) = abs( a(:,i)' * C * a(:,i) );
    BPrad(i) = abs( a1(:,i)' * R1 * a1(:,i) );
    BPcomm(i) = abs( a2(:,i)' * Wsum * a2(:,i) );
end


% Plot Beampattern
figure; plot(angleSpaceDeg,Pdesired); hold on
plot(angleSpaceDeg, Pd, 'LineWidth', 2); hold on;
plot(angleSpaceDeg, BPrad, 'LineWidth', 2); hold on;
plot(angleSpaceDeg, BPcomm, 'LineWidth', 2); hold on;
for n = 1:N
    line([commAngles(n) commAngles(n)], [0 10],'LineStyle', '--');
    hold on
end
hold off; title('Beampattern 1'); 
legend('Desired','Sep1','Rad Only','Comm Only','Comm Angles'); grid on




























