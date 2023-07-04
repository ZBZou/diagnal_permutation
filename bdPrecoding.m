function X = bdPrecoding(H,N)
% Block diagonalization (BD) precoding function
% Inputs:
%   H: channel matrix
%   P: transmit power per user
%   N: number of receive antennas
% Output:
%   X: precoded signal matrix

% Number of users
K = size(H, 2);

% Construct interference matrix
W = zeros(N, K);
for k = 1:K
    W(:, k) = null(H(:, [1:k-1, k+1:K])');
end

% Compute precoding matrix
V = zeros(K, K);
for k = 1:K
    V(:, k) = (W(:, k)'*H(:, k))/(norm(W(:, k))^2 + eps);
end

% % Scale precoding matrix by transmit power
% V = sqrt(P)*V./vecnorm(V);

% Precoded signal matrix
X = V';

end
