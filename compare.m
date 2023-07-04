%compare
clear all; 
%%
Nu=8;
M = 4;
K = log2(M);
msg_bit = randi([0 1],[K*Nu,1]); % Bit generation
%——————————— Transmitter ——————————————
symbol = qammod(msg_bit,M,'InputType', 'bit','UnitAveragePower', true);
x = symbol;
H = (randn(Nu,Nu)+1i*randn(Nu,Nu))/sqrt(2);
%%
% H = circshift(H,1);
[Q_temp,R_temp] = qr(H');
L=R_temp'; Q=Q_temp';
%% DPC
xp = x;
% pre-equlization
% xp = inv(diag(diag(L)))*x;

for m=2:Nu % GS
   xp(m,:) = xp(m,:) - L(m,1:m-1)/L(m,m)*xp(1:m-1,:);
end
Tx_signal = Q'*xp; % DPC/TH encoder
%% 2-D HOGMT
[Ur,Sr,Vr] = svd(L);
[U,S,V] = svd(H);
Tx_hogmt1 = V*pinv(S)*U'*(diag(diag(L)))*x;
Tx_hogmt2 = Q'*Vr*pinv(Sr)*Ur'*(diag(diag(L)))*x;
%%
pn= Water_Pouring(diag(S).^2,1,1);
K = sqrt(pn)*S;

