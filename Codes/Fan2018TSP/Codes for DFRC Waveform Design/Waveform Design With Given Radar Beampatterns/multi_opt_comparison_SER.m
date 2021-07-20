clc;
clear all;
close all;
warning off;
N = 16;
K = 4;
L = 20;
power = 10^(0/10);
N_montecarlo = 1000;
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
% R = waveform_mainbm_covmat( N,L,a,theta_target,theta,power );
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
    X_trdoff1 = tradeoff_comrad(rou,H,Y,power,X_orth);
    X_trdoff2 = tradeoff_comrad(rou,H,Y,power,X_arbi);
    X_trdoff3 = tradeoff_comrad_per_ant(rou,H,Y,power,X_orth);
    X_trdoff4 = tradeoff_comrad_per_ant(rou,H,Y,power,X_arbi);
%     X_trdoff5 = tradeoff_comrad_CM(rou,H,Y,power,X_arbi);
    for ii = 1:length(SNRdB)
        N0 = power/(10^(SNRdB(ii)/10));
        W = sqrt(N0/2)*(randn(K,L)+j*randn(K,L))/sqrt(2);
        Z_CR1 = H.'*X_orth+W;
        Z_CR2 = H.'*X_arbi+W;
        Z_CR3 = H.'*X_trdoff1+W;
        Z_CR4 = H.'*X_trdoff2+W;
        Z_CR5 = H.'*X_trdoff3+W;
        Z_CR6 = H.'*X_trdoff4+W;
%         Z_CR7 = H.'*X_trdoff5+W;
        Z_AWGN = Y+W;
        for mm = 1:L
            sym_demod_CR1(:,mm) = QPSK_demod(Z_CR1(:,mm));
            err_temp_CR1(mm)=sum(sym_demod_CR1(:,mm)~=Y(:,mm));
            sym_demod_CR2(:,mm) = QPSK_demod(Z_CR2(:,mm));
            err_temp_CR2(mm)=sum(sym_demod_CR2(:,mm)~=Y(:,mm));
            sym_demod_CR3(:,mm) = QPSK_demod(Z_CR3(:,mm));
            err_temp_CR3(mm)=sum(sym_demod_CR3(:,mm)~=Y(:,mm));
            sym_demod_CR4(:,mm) = QPSK_demod(Z_CR4(:,mm));
            err_temp_CR4(mm)=sum(sym_demod_CR4(:,mm)~=Y(:,mm));
            sym_demod_CR5(:,mm) = QPSK_demod(Z_CR5(:,mm));
            err_temp_CR5(mm)=sum(sym_demod_CR5(:,mm)~=Y(:,mm));
            sym_demod_CR6(:,mm) = QPSK_demod(Z_CR6(:,mm));
            err_temp_CR6(mm)=sum(sym_demod_CR6(:,mm)~=Y(:,mm));
%             sym_demod_CR7(:,mm) = QPSK_demod(Z_CR7(:,mm));
%             err_temp_CR7(mm)=sum(sym_demod_CR7(:,mm)~=Y(:,mm));
            sym_demod_AWGN(:,mm) = QPSK_demod(Z_AWGN(:,mm));
            err_temp_AWGN(mm)=sum(sym_demod_AWGN(:,mm)~=Y(:,mm));
        end
        err_CR1(ii,nn) = sum(err_temp_CR1);
        err_CR2(ii,nn) = sum(err_temp_CR2);
        err_CR3(ii,nn) = sum(err_temp_CR3);
        err_CR4(ii,nn) = sum(err_temp_CR4);
        err_CR5(ii,nn) = sum(err_temp_CR5);
        err_CR6(ii,nn) = sum(err_temp_CR6);
%         err_CR7(ii,nn) = sum(err_temp_CR7);
        err_AWGN(ii,nn) = sum(err_temp_AWGN);
    end
    clc
    disp(['Progress - ',num2str((nn-1)*length(SNRdB)+ii),'/',num2str(length(SNRdB)*N_montecarlo)]);
end
SER_CR1 = sum(err_CR1,2)/(K*L*N_montecarlo);
SER_CR2 = sum(err_CR2,2)/(K*L*N_montecarlo);
SER_CR3 = sum(err_CR3,2)/(K*L*N_montecarlo);
SER_CR4 = sum(err_CR4,2)/(K*L*N_montecarlo);
SER_CR5 = sum(err_CR5,2)/(K*L*N_montecarlo);
SER_CR6 = sum(err_CR6,2)/(K*L*N_montecarlo);
% SER_CR7 = sum(err_CR7,2)/(K*L*N_montecarlo);
SER_AWGN = sum(err_AWGN,2)/(K*L*N_montecarlo);
%%
figure(1);
semilogy(SNRdB,SER_CR1,'x-','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(SNRdB,SER_CR2,'o-','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(SNRdB,SER_CR3,'^-','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(SNRdB,SER_CR4,'*-','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(SNRdB,SER_CR5,'d--','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(SNRdB,SER_CR6,'+--','LineWidth',1.5,'MarkerSize',8);hold on;
% semilogy(SNRdB,SER_CR7,'+--','LineWidth',1.5,'MarkerSize',8);hold on;
semilogy(SNRdB,SER_AWGN,'v--','LineWidth',1.5,'MarkerSize',8);grid on;
xlabel('Transmit SNR (dB)');
ylabel('SER');
legend('Omni-Strict','Directional-Strict','Omni-Tradeoff-Total,\rho = 0.2','Directional-Tradeoff-Total,\rho = 0.2','Omni-Tradeoff-perAnt,\rho = 0.2','Directional-Tradeoff-perAnt,\rho = 0.2','Zero MUI');
figure(2);
plot(theta*180/pi,10*log10(diag(a'*X_orth*X_orth'*a)/real(trace(X_orth*X_orth'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_arbi*X_arbi'*a)/real(trace(X_arbi*X_arbi'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_trdoff1*X_trdoff1'*a)/real(trace(X_trdoff1*X_trdoff1'))),'LineWidth',1.5);hold on;
plot(theta*180/pi,10*log10(diag(a'*X_trdoff2*X_trdoff2'*a)/real(trace(X_trdoff2*X_trdoff2'))),'LineWidth',1.5);grid on;
xlim([-90,90]);
xlabel('\theta(deg)');
ylabel('Beampattern');
legend('ComRad-Orthogonal','ComRad-Arbitrary','ComRad-Tradeoff-Orth','ComRad-Tradeoff-arbi');








