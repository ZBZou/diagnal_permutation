clear all; clf
N_frame=10; N_packet=200; % Number of frames/packet and Number of packets
b=2; % Number of bits per QPSK symbol
NT=10;  NR=1;  N_user=10;  
N_pbits = N_frame*NT*b; % Number of bits in a packet
N_tbits = N_pbits*N_packet; % Number of total bits
SNRdBs = [0:2:30]; sq2=sqrt(2);
BER = zeros(length(SNRdBs),1);
for i_SNR=1:length(SNRdBs)
   SNRdB=SNRdBs(i_SNR); N_ebits=0; rand('seed',1); randn('seed',1);
   sigma2 = NT*0.5*10^(-SNRdB/10); sigma = sqrt(sigma2);
   for i_packet=1:N_packet
      msg_bit = randi([0 1],[N_pbits,1]); % bit generation
      symbol = QPSK_mapper(msg_bit).';
      x = reshape(symbol,NT,N_frame); 
      H = (randn(N_user,NT)+1i*randn(N_user,NT))/sq2;
      W = zeros(NT, N_user);
      for k = 1:N_user
          W(:, k) = null(H(:, [1:k-1, k+1:N_user])');
      end
      Tx_Data = W*x;
      Rx = H*Tx_Data + sigma*(randn(N_user,N_frame)+1i*randn(N_user,N_frame));
      W_H=H*W;
      EQ = W_H'*inv(W_H*W_H'); % Equalizer for the 1st user
    
      y = EQ*Rx;
      symbol_hat = reshape(y,NT*N_frame,1);
      symbol_sliced = QPSK_slicer(symbol_hat);
      demapped = QPSK_demapper(symbol_sliced);
      length(find(msg_bit~=demapped'))
      N_ebits = N_ebits + length(find(msg_bit~=demapped')) 
   end
   BER(i_SNR) = N_ebits/N_tbits;
end
semilogy(SNRdBs,BER,'-o'), grid on
