%%
G = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
A = [2 3 5 7; 11 13 17 19; 23 29 31 37; 41 43 47 53];
C = [1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 4];
[qa_temp,ra_temp] = qr(A');
La = ra_temp'; Qa = qa_temp';
B = G*A;
[qb_temp,rb_temp] = qr(B');
Lb = rb_temp'; Qb = qb_temp';
[qg_temp,rg_temp] = lu(G');
Lg = rg_temp'; Qg = qg_temp';
G'*C*G;
%%

% MIMO_channel_cap_ant_sel_optimal.m
clear all; clf
NT=4; NR=4; MaxIter=1000;
sel_ant=4; I=eye(NR,NR); sq2=sqrt(2);
SNRdBs=[0:2:20];
for i_SNR=1:length(SNRdBs)
   SNRdB = SNRdBs(i_SNR);
   SNR_sel_ant = 10^(SNRdB/10)/sel_ant;
   cum = 0;
   for i=1:MaxIter
      H = (randn(NR,NT)+i*randn(NR,NT))/sq2;
      if sel_ant>NT||sel_ant<1
         error('sel_ant must be between 1 and NT!');
      else
%          indices = nchoosek([1:NT],sel_ant);
         indices = perms([1:NT]);
      end
      for n=1:size(indices,1)
         Hn = H(:,indices(n,:));
         log_SH(n)=log2(real(det(I+SNR_sel_ant*Hn*Hn'))); % Eq.(12.22)
      end
      cum = cum + max(log_SH);
   end
   sel_capacity(i_SNR) = cum/MaxIter;
end
plot(SNRdBs,sel_capacity,'-ko', 'LineWidth',2); hold on;
xlabel('SNR[dB]'), ylabel('bps/Hz'), grid on;

%%
clear all;
close all;
Nu= 4;
M = 16;
K = log2(M);
T = 200;
SNRdB = 20;
sigma2 = Nu*0.5*10^(-SNRdB/10); sigma = sqrt(sigma2);
msg_bit = randi([0 1],[K*1*T,1]); % Bit generation
%——————————— Transmitter ——————————————
symbol = qammod(msg_bit,M,'InputType', 'bit','UnitAveragePower', true);
% x = reshape(symbol, Nu, T);
x = [];
for i = 1:Nu
   x(i,:) = symbol';
end
H = (randn(Nu,Nu)+1i*randn(Nu,Nu))/sqrt(2);
indices = perms([1:Nu]);
P = [];
p = [];
plt = [];
for i = 1:size(indices,1)
   Hi = H([indices(i,:)],:);
   %%
   % Hi = circshift(H,1);
   [Q_temp,R_temp] = qr(Hi');
   L=R_temp'; Q=Q_temp';
   %% DPC
   xp = x;
   % pre-equlization
   % xp = inv(diag(diag(L)))*x;

   for m=2:Nu % GS
      xp(m,:) = xp(m,:) - L(m,1:m-1)/L(m,m)*xp(1:m-1,:);
   end
   Tx_signal = Q'*xp; % DPC/TH encoder
   P(i) = abs(trace(Tx_signal'*Tx_signal)/(T*Nu));
   p(i) = max(max(Tx_signal.*conj(Tx_signal)))/P(i);
   %% 2-D HOGMT
   [Ur,Sr,Vr] = svd(L);
   [U,S,V] = svd(H);
%    Tx_hogmt1 = V*pinv(S)*U'*x;
%    Tx_hogmt2 = Q'*Vr*pinv(Sr)*Ur'*(diag(diag(L)))*x;
%    %%
%    pn= Water_Pouring(diag(S).^2,1,1);
%    K = sqrt(pn)*S;
%    Tx_hogmt1 = V*pinv(S)*U'*diag(K)*x;
   
   rx = H*Tx_signal + ...
   sigma*(randn(Nu,T)+1i*randn(Nu,T));
   y = inv(diag(diag(L)))*rx;
%    
%    Rx = H*Tx_hogmt1 + ...
%    sigma*(randn(Nu,T)+1i*randn(Nu,T));
plt(:,:,i) = Tx_signal;
plt2(:,:,i) = y;

end

[~, Pmin] = min(P);
[~, Pmax] = max(P);
[~, pmin] = min(p);
[~, pmax] = max(p);

indx = [Pmin,Pmax,pmin,pmax];
tlts = ["MinimalPower", "MximalPower","MinimalPAPR","MaximalPAPR"];

save_file = 0;

for i = 1:Nu
label = cell(1,Nu);
figure;
Tx_i = plt2(:,:,indx(i));
   for j = 1:Nu  
      scatter(real(Tx_i(j,:)),imag(Tx_i(j,:)));
      hold on;
      label{j} = "user"+num2str(j);
   end
legend(label);
xlabel('In-phase');
ylabel('Quadrature');
title(tlts(i));
filename = tlts(i);
[h, wd, ht] = tightfig();
%%

if save_file == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end
end

P = 10*log10(abs(P));
figure;
yyaxis left;
plot([1:size(indices,1)], P, '-o');
ylabel('Power [dB]');
yyaxis right;
% ylim([0 20]);
plot([1:size(indices,1)], p, '--*');
legend("AP","PAPR");
xlim([1 size(indices,1)]);
% ylim([1 4]);
xlabel('\pi(m)');
ylabel('Ratio');
filename = "Power";
if save_file == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end