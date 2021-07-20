%%Producing Fig. 3 and Fig. 4
clc;
clear all;
close all;
warning off;
N = 16;
K = 4;
L = 20;
power = 10^(0/10);
N_montecarlo = 100;
SNRdB = [-2:2:12];
%%-------------Radar Parameters-------------------
delta=pi/180;
theta=-pi/2:delta:pi/2;
theta_target=[-pi*10/180,-pi*5/180,0,pi*5/180,pi*10/180];
target_DoA=[-pi/3,0,pi/3]; 
beam_width=9;
l=ceil((target_DoA+pi/2*ones(1,length(target_DoA)))/(delta)+ones(1,length(target_DoA)));
Pd_theta=zeros(length(theta),1);
for ii=1:length(target_DoA)
    Pd_theta(l(ii)-(beam_width-1)/2:l(ii)+(beam_width-1)/2,1)=ones(beam_width,1);
end
c=3e8;
fc=3.2e9;
lamda=c/fc;
spacing=lamda/2;
for tt=1:N
    for jj=1:length(theta)
        a(tt,jj)=exp(j*pi*(tt-ceil((N)/2))*sin(theta(jj)));
    end
end

% for tt=1:N
%     for jj=1:length(theta)
%         a1(tt,jj)=exp(j*pi*(tt-ceil((N)/2))*sin(theta(jj)));
%     end
% end


R = waveform_design_multibm_covmat( Pd_theta,N,L,a,theta,power);
%R = waveform_mainbm_covmat( N,L,a,theta_target,theta,power );
F = chol(R)';
rou = 0.2;
rou1 = 0.8;
amp = sqrt(power);
for nn = 1:N_montecarlo
    H = (randn(N,K)+j*randn(N,K))/sqrt(2);
    N_pbits = 2*K*L;
    msg_bits = randi([0,1],1,N_pbits);
    Y = reshape(QPSK_mapper(msg_bits),[K,L]);
    X_orth = Orthogonal_Com_Rad( H,Y,power );
    X_arbi = Arbitrary_Com_Rad( H,Y,power,F );
    H_pinv = pinv(H.');
    tt = trace(H_pinv*Y*Y'*H_pinv');
    ff = sqrt(L*power/tt);
    X_trdoff1 = tradeoff_comrad(rou,H,Y,power,X_orth);
    X_trdoff2 = tradeoff_comrad(rou,H,Y,power,X_arbi);
    X_trdoff3 = tradeoff_comrad_per_ant(rou,H,Y,power,X_orth);
    X_trdoff4 = tradeoff_comrad_per_ant(rou,H,Y,power,X_arbi);
    for ii = 1:length(SNRdB)
        N0 = power/(10^(SNRdB(ii)/10));
        X_ZF = ff*H_pinv*Y;
        for mm = 1:L
            MUI_orth(:,mm) = abs(H.'*X_orth(:,mm)-amp*Y(:,mm)).^2;
            MUI_arbi(:,mm) = abs(H.'*X_arbi(:,mm)-amp*Y(:,mm)).^2;
            MUI_trdoff1(:,mm) = abs(H.'*X_trdoff1(:,mm)-amp*Y(:,mm)).^2;
            MUI_trdoff2(:,mm) = abs(H.'*X_trdoff2(:,mm)-amp*Y(:,mm)).^2;
            MUI_trdoff3(:,mm) = abs(H.'*X_trdoff3(:,mm)-amp*Y(:,mm)).^2;
            MUI_trdoff4(:,mm) = abs(H.'*X_trdoff4(:,mm)-amp*Y(:,mm)).^2;
            %MUI_ZF(:,mm) = abs(H.'*X_ZF(:,mm)-amp*Y(:,mm)).^2;
        end
        EMUI_orth = mean(MUI_orth,2);
        EMUI_arbi = mean(MUI_arbi,2);
        EMUI_trdoff1 = mean(MUI_trdoff1,2);
        EMUI_trdoff2 = mean(MUI_trdoff2,2);
        EMUI_trdoff3 = mean(MUI_trdoff3,2);
        EMUI_trdoff4 = mean(MUI_trdoff4,2);
        %EMUI_ZF = mean(MUI_ZF,2);
        sumrate_orth(ii,nn) = sum(log2(1+power./(EMUI_orth+N0*ones(K,1))));
        sumrate_arbi(ii,nn) = sum(log2(1+power./(EMUI_arbi+N0*ones(K,1))));
        sumrate_trdoff1(ii,nn) = sum(log2(1+power./(EMUI_trdoff1+N0*ones(K,1))));
        sumrate_trdoff2(ii,nn) = sum(log2(1+power./(EMUI_trdoff2+N0*ones(K,1))));
        sumrate_trdoff3(ii,nn) = sum(log2(1+power./(EMUI_trdoff3+N0*ones(K,1))));
        sumrate_trdoff4(ii,nn) = sum(log2(1+power./(EMUI_trdoff4+N0*ones(K,1))));
        %sumrate_ZF(ii,nn) = sum(log2(1+power./(EMUI_ZF+N0*ones(K,1))));
        sumrate_lim(ii,nn) = sum(log2(1+power./(N0*ones(K,1))));
    end
    clc
    disp(['Progress - ',num2str((nn-1)*length(SNRdB)+ii),'/',num2str(length(SNRdB)*N_montecarlo)]);
end
%%
figure(1);
plot(SNRdB,mean(sumrate_orth,2),'x-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_arbi,2),'o-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_trdoff1,2),'^-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_trdoff2,2),'*-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_trdoff3,2),'d--','LineWidth',1.5,'MarkerSize',8);hold on;
plot(SNRdB,mean(sumrate_trdoff4,2),'+--','LineWidth',1.5,'MarkerSize',8);hold on;
%plot(SNRdB,mean(sumrate_ZF,2),'s-','LineWidth',1.5,'MarkerSize',8);
plot(SNRdB,mean(sumrate_lim,2),'v--','LineWidth',1.5,'MarkerSize',8);
grid on;
xlabel('Transmit SNR (dB)');
ylabel('Average Achievable Sum Rate (bps/Hz)');
legend('Omni-Strict','Directional-Strict','Omni-Tradeoff-Total,\rho = 0.2','Directional-Tradeoff-Total,\rho = 0.2','Omni-Tradeoff-perAnt,\rho = 0.2','Directional-Tradeoff-perAnt,\rho = 0.2','Zero MUI');
figure(2);
plot(theta*180/pi,10*log10(diag(a'*X_orth*X_orth'*a)/real(trace(X_orth*X_orth'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_arbi*X_arbi'*a)/real(trace(X_arbi*X_arbi'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_trdoff1*X_trdoff1'*a)/real(trace(X_trdoff1*X_trdoff1'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_trdoff2*X_trdoff2'*a)/real(trace(X_trdoff2*X_trdoff2'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_trdoff3*X_trdoff3'*a)/real(trace(X_trdoff3*X_trdoff3'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_trdoff4*X_trdoff4'*a)/real(trace(X_trdoff4*X_trdoff4'))),'LineWidth',1.5);grid on;
xlim([-90,90]);
xlabel('\theta(deg)');
ylabel('Beampattern');
legend('Omni-Strict','Directional-Strict','Omni-Tradeoff-Total,\rho = 0.2','Directional-Tradeoff-Total,\rho = 0.2','Omni-Tradeoff-perAnt,\rho = 0.2','Directional-Tradeoff-perAnt,\rho = 0.2');








