% Simulating the DFRC Beamformer presented in 
%   "Joint Radar and Communication Design: Applications, State-of-the-art,
%    and the Road Ahead" - Fan Liu et. al. 2019

% Nate Raymondi, 1/23/2020

% Assign frequencies and propagation speed
Nt = 64;             % number of transmit antennas (BS)
Nr = 10;             % number of receiver antennas (UE)
fc = 24e9;
c = 3e8;
lamb = c/fc;
K = 8; 
L = 4;
P = 1;

% create antenna position arrays, both lambda/2 spacing
spacing = lamb/2;
TxAntLoc = spacing*[0:Nt-1];
RxAntLoc = spacing*[0:Nr-1];

% set up angles, assume L = 4, K = 8
%   AoDs (from BS to L scatterers)
AoD1 = [-30;0]; AoD2 = [-10;0]; AoD3 = [8;0]; AoD4 = [25;0];
%   AoAs (from L scatterers to UE)
AoA1 = [-30;0]; AoA2 = [-10;0]; AoA3 = [8;0]; AoA4 = [25;0];
%   Radar Target reflections
tgt1 = [-75;0]; tgt2 = [-45;0]; tgt3 = [52;0]; tgt4 = [80;0];

Kangles = [tgt1, tgt2, AoD1, AoD2, AoD3, AoD4, tgt3, tgt4];
targets = [tgt1, tgt2, tgt3, tgt4];
Langles = [AoD1, AoD2, AoD3, AoD4];
LanglesUE = [AoA1, AoA2, AoA3, AoA4];

% define steering vectors
TxArray = phased.ULA('NumElements',Nt, 'ElementSpacing', spacing);
TxSteerVec = phased.SteeringVector('SensorArray',TxArray);
RxArray = phased.ULA('NumElements',Nr, 'ElementSpacing', spacing);
RxSteerVec = phased.SteeringVector('SensorArray',RxArray);

% a(\theta) - all K targets
A_K = zeros(Nt,K); 
for k = 1:K
    A_K(:,k) = TxSteerVec(fc,Kangles(:,k));
end

% a(\psi) - only AoDs
A_L = zeros(Nt,L); B_L = zeros(Nr,L);
for i = 1:L
    A_L(:,i) = TxSteerVec(fc,Langles(:,i));
    B_L(:,i) = RxSteerVec(fc,LanglesUE(:,i));
end

% draw \beta_l from standard complex gaussian distribution (a+jb), a,b~N(0,1/2)
beta = randn([L,1]) + 1i*randn([L,1]);
betaMat = diag(beta);

% Create the channel matrix according to Eq. (5)
term = zeros(L,Nr,Nt);
for i = 1:L
    term(i,:,:) = beta(i)*B_L(:,i)*transpose(A_L(:,i));
end
H = squeeze(sum(term));

% channel approximation \tilde{H} = diag(beta)A'(theta_1)
H_approx = betaMat*transpose(A_L);

% create ZF beamformers Eq. (38)
F_bs = H_approx' * inv(H_approx * H_approx');
W_ue = inv(B_L' * B_L) * B_L';

% Analog Beamformer F_rf Eq. (40)
F_rf = zeros(Nt,K);
for i = 1:K
    F_rf(:,i) = conj(A_K(:,i));
end

% Auxilary Matrix F_aux, where H_approx * F_aux = 0
nullMat = null(H_approx);
F_aux_NSP = nullMat(:,1:K-L);

% Constructing the Digital Beamformer F_bb Eq. (42,43)
svdMat = F_rf' * [F_bs,F_aux_NSP];
[U,S,V] = svd(svdMat);
F_bb = sqrt(P/(K*Nt)) * U * V';

% Covariance Matrix Eq. (39)
R = F_rf * F_bb * F_bb' * F_rf';    R = normalize(R,'norm');
R_zf = F_bs * F_bs';                R_zf = normalize(R_zf,'norm');

% Plot the beampatterns d(\theta) = transpose(a(\theta)) * R * conj(a(\theta))
% Find ALL steering vectors on [-pi,pi] to plot the beampattern 
angleSpace = linspace(-pi/2, pi/2, 360);
angleSpaceDeg = linspace(-90, 90, 360);
a = zeros(Nt, length(angleSpace));
for j = 1:Nt
    a(j,:) = exp((1i * 2 * pi * TxAntLoc(j) / lamb) .* sin(angleSpace));
end

 dZF = zeros(size(angleSpace)); d = zeros(size(angleSpace));
for i = 1:length(angleSpace)
    dZF(i) = transpose(a(:,i)) * R_zf * conj(a(:,i));
    d(i) = transpose(a(:,i)) * R * conj(a(:,i));
end


% PLOTTING BEAMPATTERNS -----------------------------------------------------
tgts_plot = targets(1,:);           % quick fix to deal with some annoying 
L_angles_plot = Langles(1,:);       % matrix dimensions

figure
p = plot(angleSpaceDeg, mag2db(abs(dZF)), 'LineWidth', 2);
hold on
for i = 1:length(tgts_plot)
    p_tgt = line([tgts_plot(i) tgts_plot(i)], [min(mag2db(abs(dZF)))...
        max(mag2db(abs(dZF)))], 'Color', 'black', 'LineStyle', '--');
    hold on
end
hold on
for i = 1:length(L_angles_plot)
    p_c = line([L_angles_plot(i) L_angles_plot(i)], [min(mag2db(abs(dZF)))...
        max(mag2db(abs(dZF)))], 'Color', 'magenta', 'LineStyle', '--');
    hold on
end
xlabel('Angle Space [-90^\circ,90^\circ]'); ylabel('Magnitude (dB)')
title('BS ZF Transmit Beampattern'); grid on; axis tight
set(gcf,'color','w'); set(gcf, 'Position',  [50, 100, 1000, 400])
legend([p,p_tgt,p_c],'ZF Beampattern','Radar Directions','Comms Directions');

figure
p = plot(angleSpaceDeg, mag2db(abs(d)), 'LineWidth', 2);
hold on
for i = 1:length(tgts_plot)
    p_tgt = line([tgts_plot(i) tgts_plot(i)], [min(mag2db(abs(d)))...
        max(mag2db(abs(d)))], 'Color', 'black', 'LineStyle', '--');
    hold on
end
hold on
for i = 1:length(L_angles_plot)
    p_c = line([L_angles_plot(i) L_angles_plot(i)], [min(mag2db(abs(d)))...
        max(mag2db(abs(d)))], 'Color', 'magenta', 'LineStyle', '--');
    hold on
end
xlabel('Angle Space [-90^\circ,90^\circ]'); ylabel('Magnitude (dB)')
title('BS DFRC Transmit Beampattern'); grid on; axis tight
set(gcf,'color','w'); set(gcf, 'Position',  [50, 100, 1000, 400])
legend([p,p_tgt,p_c],'DFRC Beampattern','Radar Directions','Comms Directions');


