%%Producing Fig. 5
clc;
clear all;
close all;
warning off;
N = 16;
% K = 4;
L = 20;
power = 10^(0/10);
amp = sqrt(power);
N_montecarlo = 100;
SNRdB = 10;
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
SNRr = 10^(-6/10);
uu = 36;

% H = (randn(N,K)+j*randn(N,K))/sqrt(2);
% N_pbits = 2*K*L;
% msg_bits = randint(1,N_pbits);
% Y = reshape(QPSK_mapper(msg_bits),[K,L]);
% X1 = sqrt(N)*Orthogonal_Com_Rad( H,Y,power );
% RMSE = CRB_Orthogonal( X1,a(:,uu),theta(uu),SNR );

Nii = 20;
N0 = power/(10^(SNRdB/10));
Nkk  = 3;
for kk = 1:Nkk
    K= 4+(kk-1)*2;
    for nn = 1:N_montecarlo
        H = (randn(N,K)+j*randn(N,K))/sqrt(2);
        N_pbits = 2*K*L;
        msg_bits = randi([0,1],1,N_pbits);
        Y = reshape(QPSK_mapper(msg_bits),[K,L]);
        X1 = Orthogonal_Com_Rad( H,Y,power );
%         RMSE1 = CRB_Orthogonal( X1,a(:,uu),theta(uu),SNRr );
%         H_pinv = pinv(H.');
%         tt = trace(H_pinv*Y*Y'*H_pinv');
%         X3 = sqrt(N*power/tt)*H_pinv*Y;
        for ii = 1:Nii-1
            rou = ii/Nii;
            X2 = sqrt(N)*tradeoff_comrad(rou,H,Y,power,X1);
            %         for mm = 1:L
            %             MUI1(:,mm) = abs(H.'*X1(:,mm)/sqrt(N)-amp*Y(:,mm)).^2;
            %             MUI2(:,mm) = abs(H.'*X2(:,mm)/sqrt(N)-amp*Y(:,mm)).^2;
            %             MUI3(:,mm) = abs(H.'*X3(:,mm)/sqrt(N)-amp*Y(:,mm)).^2;
            %         end
%             MUI1 = abs(H.'*X1/sqrt(N)-amp*Y).^2;
            MUI2 = abs(H.'*X2/sqrt(N)-amp*Y).^2;
%             MUI3 = abs(H.'*X3/sqrt(N)-amp*Y).^2;
%             EMUI1 = mean(MUI1,2);
            EMUI2 = mean(MUI2,2);
%             EMUI3 = mean(MUI3,2);
%             sumrate1(ii,kk,nn) = sum(log2(1+power./(EMUI1+N0*ones(K,1))));
            sumrate2(ii,kk,nn) = sum(log2(1+power./(EMUI2+N0*ones(K,1))))/K;
%             sumrate3(ii,kk,nn) = sum(log2(1+power./(EMUI3+N0*ones(K,1))));
%             RMSE1(ii,kk,nn) = CRB_Orthogonal( X1,a(:,uu),theta(uu),SNRr );
            PD2(ii,kk,nn) = PD_Orthogonal( X2,a(:,uu),SNRr );
%             RMSE3(ii,kk,nn) = CRB_Orthogonal( X3,a(:,uu),theta(uu),SNRr );
            clc
            disp(['Progress - ',num2str((kk-1)*N_montecarlo*Nii+(nn-1)*Nii+ii),'/',num2str(Nii*N_montecarlo*Nkk)]);
        end
    end
end
%%
figure(1);
% plot(mean(sumrate2,2),mean(RMSE1,2),'x-','LineWidth',1.5,'MarkerSize',8);hold on;
for kk = 1:Nkk
    plot(mean(sumrate2(:,kk,:),3),mean(PD2(:,kk,:),3),'-','LineWidth',1.5,'MarkerSize',8);hold on;
end
% plot(mean(sumrate2,2),mean(RMSE3,2),'^-','LineWidth',1.5,'MarkerSize',8);hold on;
grid on;
legend('K = 4','K = 6','K = 8');
xlabel('Average achievable rate (bps/Hz/user)');
ylabel('P_D');



