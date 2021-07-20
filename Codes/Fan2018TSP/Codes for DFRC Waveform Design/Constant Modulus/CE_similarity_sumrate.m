%% Producing Fig. 9 and Fig. 10
clc;
clear all;
close all;
warning off;
N = 16;
K = 4;
N_montecarlo = 50;
power = 10^(0/10);
amp = sqrt(power/N);
SNRdB = 10;

cle = 1;

%%-------------Radar Parameters-------------------
fc = 10;
B = 20;
Tp = 2;
X0 = chirp_generator( N,fc,Tp,B );
L = length(X0(1,:));


N_pbits = 2*K*L;
Nii = 11;
%%
N0 = power/(10^(SNRdB/10));
for jj = 1:N_montecarlo
    H = (randn(N,K)+j*randn(N,K))/sqrt(2);
    amp = sqrt(power/N);
    H_wave = amp*[real(H),imag(H);-imag(H),real(H)]; % Equivalent Real Channel
    for ii = 1:Nii
        ee = 0.05+(ii-1)*0.19;
%         ee = 0.4;
        delta = acos(1-ee^2/2);
        msg_bits = randi([0,1],1,N_pbits);
        Y = reshape(QPSK_mapper(msg_bits),[K,L]);
        try
            parfor nn = 1:L
                warning off;
                y = Y(:,nn);
                x0 = X0(:,nn);
                y_wave = sqrt(power*cle)*[real(y);imag(y)];
                %         noise = N0*(randn(K,1)+j*randn(K,1))/sqrt(2);
                l = angle(x0)-delta*ones(N,1);
                u = angle(x0)+delta*ones(N,1);     %Upper and Lower Bound
                [xr_lb,~] = QCQP_LB1( H_wave,y_wave,N,l,u);
                x_lb = xr_lb(1:N)+j*xr_lb(N+1:2*N);
                X_lb(:,nn) = amp*x_lb;
                MUI_lb(:,nn) = abs(amp*H.'*x_lb - sqrt(power*cle)*y).^2;
                %         [xr_nb,nbd(ii,nn,jj)] = normalize_UB( H_wave,y_wave,xr_lb,N,l,u);
                %         x_nb = xr_nb(1:N)+j*xr_nb(N+1:2*N);
                %         X_nb(:,nn) = amp*x_nb;
                %         MUI_nb(:,nn) = abs(amp*H.'*x_nb - sqrt(power*cle)*y).^2;
                %         [xr_ub,ubd(ii,nn,jj)] = QCQP_UB( H_wave,y_wave,N,l,u,xr_nb); % fmincon UB
                %         x_ub = xr_ub(1:N)+j*xr_ub(N+1:2*N);
                %         X_ub(:,nn) = amp*x_ub;
                %         MUI_ub(:,nn) = abs(amp*H.'*x_ub - sqrt(power*cle)*y).^2;
                [x_opt,~] = CE_similarity_ComRad( H,y,power,ee,x0,cle);
                X_opt(:,nn) = amp*x_opt;
                MUI_opt(:,nn) = abs(amp*H.'*x_opt - sqrt(power*cle)*y).^2;
%                 [x_bm,~] = CE_similarity_ComRad_benchmark( H,y,power,ee,x0,cle );
                        try
                            [x_bm,~] = CE_similarity_ComRad_benchmark( H,y,power,ee,x0,cle );
                        catch
                            [x_bm1,~] = normalize_UB( H_wave,y_wave,xr_lb,N,l,u);
                            x_bm = x_bm1(1:N)+j*x_bm1(N+1:2*N);
                        end
                X_bm(:,nn) = amp*x_bm;
                MUI_bm(:,nn) = abs(amp*H.'*x_bm - sqrt(power*cle)*y).^2;
                clc
                disp(['Progress - ',num2str((jj-1)*Nii*L+(ii-1)*L+nn),'/',num2str(L*N_montecarlo*Nii)]);
            end
        catch
            for nn = 1:L
                warning off;
                y = Y(:,nn);
                x0 = X0(:,nn);
                y_wave = sqrt(power*cle)*[real(y);imag(y)];
                l = angle(x0)-delta*ones(N,1);
                u = angle(x0)+delta*ones(N,1);     %Upper and Lower Bound
                [xr_lb,~] = QCQP_LB1( H_wave,y_wave,N,l,u);
                x_lb = xr_lb(1:N)+j*xr_lb(N+1:2*N);
                X_lb(:,nn) = amp*x_lb;
                MUI_lb(:,nn) = abs(amp*H.'*x_lb - sqrt(power*cle)*y).^2;
                [x_opt,~] = CE_similarity_ComRad( H,y,power,ee,x0,cle);
                X_opt(:,nn) = amp*x_opt;
                MUI_opt(:,nn) = abs(amp*H.'*x_opt - sqrt(power*cle)*y).^2;
                %         [x_bm,bm(ii,nn,jj)] = CE_similarity_ComRad_benchmark( H,y,power,ee,x0,cle );
                try
                    [x_bm,~] = CE_similarity_ComRad_benchmark( H,y,power,ee,x0,cle );
                catch
                    [x_bm1,~] = normalize_UB( H_wave,y_wave,xr_lb,N,l,u);
                    x_bm = x_bm1(1:N)+j*x_bm1(N+1:2*N);
                end
                X_bm(:,nn) = amp*x_bm;
                MUI_bm(:,nn) = abs(amp*H.'*x_bm - sqrt(power*cle)*y).^2;
                clc
                disp(['Progress - ',num2str((jj-1)*Nii*L+(ii-1)*L+nn),'/',num2str(L*N_montecarlo*Nii)]);
            end
        end
        
        
        C_opt(ii,jj) = K*log2(1+power*cle/N0);
        EMUI_lb = mean(MUI_lb,2);
%         EMUI_nb = mean(MUI_nb,2);
%         EMUI_ub = mean(MUI_ub,2);
        EMUI_opt = mean(MUI_opt,2);
        EMUI_bm = mean(MUI_bm,2);
        sumrate_lb(ii,jj) = sum(log2(1+power*cle./(EMUI_lb+N0*ones(K,1))));
%         sumrate_nb(ii,jj) = sum(log2(1+power*cle./(EMUI_nb+N0*ones(K,1))));
%         sumrate_ub(ii,jj) = sum(log2(1+power*cle./(EMUI_ub+N0*ones(K,1))));
        sumrate_opt(ii,jj) = sum(log2(1+power*cle./(EMUI_opt+N0*ones(K,1))));
        sumrate_bm(ii,jj) = sum(log2(1+power*cle./(EMUI_bm+N0*ones(K,1))));
        PSLR_lb(ii,jj) = 10^(PSLRindB(X_lb(1,:)/amp)/20);
%         PSLR_nb(ii,jj) = 10^(PSLRindB(X_nb(1,:)/amp)/20);
%         PSLR_ub(ii,jj) = 10^(PSLRindB(X_ub(1,:)/amp)/20);
        PSLR_opt(ii,jj) = 10^(PSLRindB(X_opt(1,:)/amp)/20);
        PSLR_bm(ii,jj) = 10^(PSLRindB(X_bm(1,:)/amp)/20);
        sf_opt(:,ii,jj) = 10.^(pc_freqwin(X_opt(1,:)/amp)/20).';
        sf_bm(:,ii,jj) = 10.^(pc_freqwin(X_bm(1,:)/amp)/20).';
    end
end

%%
% MSE_lb = mean(mean(lbd,3),2);
% MSE_nb = mean(mean(nbd,3),2);
% MSE_ub = mean(mean(ubd,3),2);
% MSE_opt = mean(mean(opt,3),2);
% MSE_bm = mean(mean(bm,3),2);
EPSLR_lb = 20*log10(mean(PSLR_lb(:,1:jj),2));
% EPSLR_nb = 20*log10(mean(PSLR_nb(:,1:jj),2));
% EPSLR_ub = 20*log10(mean(PSLR_ub(:,1:jj),2));
EPSLR_opt = 20*log10(mean(PSLR_opt(:,1:jj),2));
EPSLR_bm = 20*log10(mean(PSLR_bm(:,1:jj),2));
%%
% x0 = X0(1,:);
% x1 = X_opt(1,:)/amp;
% PSLR0 = PSLRindB( x0 );
% PSLR1 = PSLRindB( x1 );
% f0 = pc_freqwin( x0 );
% f1 = pc_freqwin( x1 );


%%
figure(1);
plot(0.2*[0:Nii-1],mean(sumrate_lb(:,1:jj-2),2),'x-','LineWidth',1.5,'MarkerSize',8);hold on;
% plot(0.2*[1:Nii],mean(sumrate_nb(:,1:jj),2),'*-','LineWidth',1.5,'MarkerSize',8);hold on;
% plot(0.2*[1:Nii],mean(sumrate_ub(:,1:jj),2),'o-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(0.2*[0:Nii-1],mean(sumrate_opt(:,1:jj-2),2),'^-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(0.2*[0:Nii-1],mean(sumrate_bm(:,1:jj-2),2),'s-','LineWidth',1.5,'MarkerSize',8);hold on;
plot(0.2*[0:Nii-1],mean(C_opt(:,1:jj-2),2),'o--','LineWidth',1.5,'MarkerSize',8);grid on;
% legend('QCQP-LB','Normalized-LB','QCQP-UB','Global Optimum','SQR-BS (benchmark)','Capacity UB');
legend('QCQP Convex Bound','Global Optimum (Proposed)','SQR-BS (Benchmark)','AWGN Capacity');
xlabel('\epsilon');
ylabel('Average achievable sum rate (bps/Hz)');
%%
% figure(2);
% semilogy(0.2*[1:Nii],MSE_lb,'x-','LineWidth',1.5);hold on;
% % semilogy(0.2*[1:Nii],MSE_nb,'*-','LineWidth',1.5);hold on;
% % semilogy(0.2*[1:Nii],MSE_ub,'o-','LineWidth',1.5);hold on;
% semilogy(0.2*[1:Nii],MSE_opt,'^-','LineWidth',1.5);hold on;
% semilogy(0.2*[1:Nii],MSE_bm,'^-','LineWidth',1.5);grid on;
% legend('QCQP-LB','Global Optimum','SQR-BS (benchmark)');
% xlabel('\epsilon');
% ylabel('MSE');

figure(3);
plot(EPSLR_lb,mean(sumrate_lb(:,1:jj),2),'x-','LineWidth',1.5);hold on;
% plot(EPSLR_nb,mean(sumrate_nb(:,1:jj),2),'*-','LineWidth',1.5);hold on;
% plot(EPSLR_ub,mean(sumrate_ub(:,1:jj),2),'o-','LineWidth',1.5);hold on;
plot(EPSLR_opt,mean(sumrate_opt(:,1:jj),2),'^-','LineWidth',1.5);hold on;
plot(EPSLR_opt,mean(sumrate_bm(:,1:jj),2),'^-','LineWidth',1.5);grid on;
legend('QCQP-LB','Global Optimum','SQR-BS (benchmark)');
xlabel('Average PSLR(dB)');
ylabel('Average achievable sum rate (bps/Hz)');
%%
sf_orig = pc_freqwin(X0(1,:)).';
figure(4);
plot(20*log10(mean(sf_opt(:,3,1:jj-1),3)),'LineWidth',1.5);hold on;
plot(20*log10(mean(sf_bm(:,3,1:jj-1),3)),'LineWidth',1.5);hold on;
plot(sf_orig,'LineWidth',1.5);grid on;
legend('Global Optimum (Proposed), \epsilon = 0.2','SQR-BS (Benchmark), \epsilon = 0.2','Chirp');
xlabel('IFFT Index');
ylabel('Pulse Comperession Gain (dB)');







