% dpc, hogmt...
clear all; 
close all;
mode=0; % Set 0/1 for Dirty or TH precoding
N_frame=1; N_packet=1000; % Number of frames/packet and Number of packets
% b=2; % Number of bits per QPSK symbol
mod = [4, 16, 64, 128];
NT=10; N_user=10; N_act_user=10; I=eye(N_act_user,NT);

SNRdBs=[0:30]; sq2=sqrt(2);
BERd = zeros(length(mod), length(SNRdBs));
BERh = zeros(length(mod), length(SNRdBs));
BERt = zeros(length(mod), length(SNRdBs));
BERz = zeros(length(mod), length(SNRdBs));
BERm = zeros(length(mod), length(SNRdBs));
BERb = zeros(length(mod), length(SNRdBs));
for k = 1:length(mod)
   M_qam = mod(k);
   b = log2(M_qam);
   N_pbits = N_frame*NT*b; % Number of bits in a packet
   N_tbits = N_pbits*N_packet; % Number of total bits
   for i_SNR=1:length(SNRdBs)
      SNRdB = SNRdBs(i_SNR); 
      N_ebitsd = 0; N_ebitst = 0;N_ebitsh = 0;N_ebitsz = 0;
      N_ebitsm = 0;N_ebitsb = 0;
      rand('seed',1); randn('seed',1);
      sigma2 = NT*0.5*10^(-SNRdB/10); sigma = sqrt(sigma2);
      for i_packet=1:N_packet
%          N_ebits1
         msg_bit = randi([0 1],[N_pbits,1]); % Bit generation
         %——————————— Transmitter ——————————————
         symbol = qammod(msg_bit,M_qam,'InputType', 'bit','UnitAveragePower', true);
         x = symbol;
         H = (randn(N_user,NT)+1i*randn(N_user,NT))/sq2;
         Combinations = nchoosek([1:N_user],N_act_user)'; % User selection
         for i=1:size(Combinations,2)
            H_used = H(Combinations(:,i),:);
            [Q_temp,R_temp] = qr(H_used);
            % Diagonal entries of R_temp are real
            minimum_l(i) = min(diag(R_temp));
         end
         [max_min_l,Index] = max(minimum_l);
         H_used = H(Combinations(:,Index),:);
         [Q_temp,R_temp] = qr(H_used');
         L=R_temp'; Q=Q_temp';
         %————————— precoding ——————————
         temp_Wz = H_used'*inv(H_used*H_used');
         temp_Wm = H_used'*inv(H_used*H_used'+sigma2*I);
         %BD
            Wb = zeros(NT, N_user);
         for nu = 1:N_user
            Wb(:, nu) = null(H_used(:, [1:nu-1, nu+1:N_user])');
         end
         xp_b = Wb*x;
      
         xp_z = temp_Wz*x;
         xp_m = temp_Wm*x;
         % pre-equlization
         xp_d = inv(diag(diag(L)))*x;
         xp_t = inv(diag(diag(L)))*x;
         % Dirty precoding
         for m=2:NT % Eqs.(13.39)(13.41)
            xp_d(m,:) = xp_d(m,:) - L(m,1:m-1)/L(m,m)*xp_d(1:m-1,:);
         end
         % TH precoding
         for m=2:NT % Eqs.(13.52)(13.54)
         xp_t(m,:) = modulo(xp_t(m,:)-L(m,1:m-1)/L(m,m)*xp_t(1:m-1,:),sq2);
         end
         
         Tx_signal_d = Q'*xp_d; % DPC
         Tx_signal_t = Q'*xp_t; 
         Tx_signal_z = xp_z; 
         Tx_signal_m = xp_m; 
         Tx_signal_b = xp_b; 
         % 2-D HOGMT
         [U,temp_sig,V] = svd(H_used);
         Tx_hogmt = V*pinv(temp_sig)*U'*x;
   %       Tx_hogmt = V*pinv(temp_sig)*U'*diag(diag(L))*x;
         %—————————— Channel and Noise ——————————————
         Rx_signal_d = H_used*Tx_signal_d + ...
         sigma*(randn(N_act_user,N_frame)+1i*randn(N_act_user,N_frame));
         Rx_signal_t = H_used*Tx_signal_t + ...
         sigma*(randn(N_act_user,N_frame)+1i*randn(N_act_user,N_frame));
         Rx_signal_z = H_used*Tx_signal_z + ...
         sigma*(randn(N_act_user,N_frame)+1i*randn(N_act_user,N_frame));
         Rx_signal_m = H_used*Tx_signal_m + ...
         sigma*(randn(N_act_user,N_frame)+1i*randn(N_act_user,N_frame));
         Rx_signal_b = H_used*Tx_signal_b + ...
         sigma*(randn(N_act_user,N_frame)+1i*randn(N_act_user,N_frame));
      
         Rx_signal_h = H_used*Tx_hogmt + ...
         sigma*(randn(N_act_user,N_frame)+1i*randn(N_act_user,N_frame));
         %—————————— Receiver ——————————————
         %y = inv(diag(diag(L)))*Rx_signal;
         yd = Rx_signal_d;
         yt = Rx_signal_t;
         yh = Rx_signal_h;
         yz = Rx_signal_z;
         ym = Rx_signal_m;
         W_H=H*Wb;
         EQ = W_H'*inv(W_H*W_H');
         yb = EQ*Rx_signal_b;
         symbol_hatd = reshape(yd,NT*N_frame,1);
         symbol_hatz = reshape(yz,NT*N_frame,1);
         symbol_hatm = reshape(ym,NT*N_frame,1);
         symbol_hath = reshape(yh,NT*N_frame,1);
         symbol_hatt = reshape(yt,NT*N_frame,1);
         symbol_hatt = modulo(symbol_hatt,sq2);
         symbol_hatb = reshape(yb,NT*N_frame,1);
         
         demapped_b = qamdemod(symbol_hatb,M_qam,'OutputType', 'bit', 'UnitAveragePower', true);
         N_ebitsb = N_ebitsb + length(find(msg_bit~=demapped_b));
         
         
         demapped_z = qamdemod(symbol_hatz,M_qam,'OutputType', 'bit', 'UnitAveragePower', true);
         N_ebitsz = N_ebitsz + length(find(msg_bit~=demapped_z));
         
         demapped_m = qamdemod(symbol_hatm,M_qam,'OutputType', 'bit', 'UnitAveragePower', true);
         N_ebitsm = N_ebitsm + length(find(msg_bit~=demapped_m));
         
         demapped_d = qamdemod(symbol_hatd,M_qam,'OutputType', 'bit', 'UnitAveragePower', true);
         N_ebitsd = N_ebitsd + length(find(msg_bit~=demapped_d));
         demapped_t = qamdemod(symbol_hatt,M_qam,'OutputType', 'bit', 'UnitAveragePower', true);
         N_ebitst = N_ebitst + length(find(msg_bit~=demapped_t));
         demapped_h = qamdemod(symbol_hath,M_qam,'OutputType', 'bit', 'UnitAveragePower', true);
         N_ebitsh = N_ebitsh + length(find(msg_bit~=demapped_h));
      end
      BERb(k,i_SNR) = N_ebitsb/N_tbits;
      BERz(k,i_SNR) = N_ebitsz/N_tbits;
      BERm(k,i_SNR) = N_ebitsm/N_tbits;
      BERd(k,i_SNR) = N_ebitsd/N_tbits;
      BERt(k,i_SNR) = N_ebitst/N_tbits;
      BERh(k,i_SNR) = N_ebitsh/N_tbits;
   end
end
%%
figure();
mesh(abs(H_used))
%%
% figure();
% semilogy(SNRdBs,BER1(1,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(1,:),'--'), grid on
% hold on;
% semilogy(SNRdBs,BER1(2,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(2,:),'--'), grid on
% hold on;
% semilogy(SNRdBs,BER1(3,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(3,:),'--'), grid on
% hold on;
% semilogy(SNRdBs,BER1(4,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(4,:),'--'), grid on
% hold on;
% xlim([0 20]);
% legend("Conventional DPC with QPSK","Proposed method with QPSK",...
%    "Conventional DPC with 16-QAM","Proposed method with 16-QAM",...
%    "Conventional DPC with 64-QAM","Proposed method with 64-QAM",...
%    "Conventional DPC with 128-QAM","Proposed method with 128-QAM",'location','southwest')
% ylabel("BER [dB]");
% xlabel("SNR [dB]")
%%
figure();
semilogy(SNRdBs,BERh(1,:),'-o'), grid on
hold on;
semilogy(SNRdBs,BERh(2,:),'-o'), grid on
hold on;
semilogy(SNRdBs,BERh(3,:),'-o'), grid on
hold on;
semilogy(SNRdBs,BERh(4,:),'-o'), grid on
hold on;
xlim([0 20]);
legend("Proposed method with QPSK",...
   "Proposed method with 16-QAM",...
   "Proposed method with 64-QAM",...
   "Proposed method with 128-QAM",'location','southwest')
ylabel("BER [dB]");
xlabel("SNR [dB]")
%%
figure();
for i = 1
   semilogy(SNRdBs,BERb(i,:),'-o'), grid on
   hold on;
   semilogy(SNRdBs,BERt(i,:),'-*'), grid on
   hold on;
%    semilogy(SNRdBs,BERz(i,:),'-o'), grid on
%    hold on;
   semilogy(SNRdBs,BERm(i,:),'-o'), grid on
   hold on;
   semilogy(SNRdBs,BERd(i,:),'-o'), grid on
   hold on;
   semilogy(SNRdBs,BERh(i,:),'-s'), grid on
   hold on;

end
% semilogy(SNRdBs,BER1(2,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(2,:),'--'), grid on
% hold on;
% semilogy(SNRdBs,BER1(3,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(3,:),'--'), grid on
% hold on;
% semilogy(SNRdBs,BER1(4,:),'-o'), grid on
% hold on;
% semilogy(SNRdBs,BER2(4,:),'--'), grid on
% hold on;
xlim([0 20]);
legend("BD with QPSK","THP with QPSK","MMSE with QPSK","DPC with QPSK","Proposed method with QPSK",...
   'location','southwest')

% legend("Conventional DPC with QPSK","Proposed method with QPSK",...
%    "Conventional DPC with 16-QAM","Proposed method with 16-QAM",...
%    "Conventional DPC with 64-QAM","Proposed method with 64-QAM",...
%    "Conventional DPC with 128-QAM","Proposed method with 128-QAM",'location','southwest')
ylabel("BER [dB]");
xlabel("SNR [dB]")
%%
set(gca, 'fontsize', 20)
filename = "BER_multi";
[h, wd, ht] = tightfig();
save_file = 1;

if save_file == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end