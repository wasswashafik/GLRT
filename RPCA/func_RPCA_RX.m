function [ RPCA_RX_out ] = func_RPCA_RX( hsi,lambda )

[N,NN,M]=size(hsi);
hsi = reshape(hsi,N*NN,M);


% max_y = max(max(hsi));
% min_x = min(min(hsi));
% hsi = (hsi-min_x)/(max_y-min_x);
% 
%% Each band is normalized 0-1
% hsi = scale_new(hsi);  % column normalization


% RPCA_MC m x n matrix of observations/data£¬Behavior dimension or number of samples, detection results are the same
% hsi = hsi';
% [A_low E_spare ] = inexact_alm_rpca_mc(hsi, lambda);    % Too slow
% E_spare = E_spare';

% % ALM Method£¬m x n matrix of observations/data  n is the smallest, that is, listed as dimension
[A_low,E_spare,iter] = inexact_alm_rpca(hsi, lambda);

%% Post-detection algorithm: RX detection on sparse matrix (classic post-processing method£ºRPCA_RX)
cov_matrix = cov(E_spare,1);  % divided by N£¨N-1£©
inv_cov=inv(cov_matrix);
xmean = mean(E_spare,1);
% E_spare_img = reshape(E_spare,N,NN,M);
RPCA_RX_out = zeros(N*NN,1);

for i = 1:N*NN
    value_test = E_spare(i,:);
    RPCA_RX_out(i,:)=(value_test-xmean)*inv_cov*(value_test-xmean)';
end
RPCA_RX_out = reshape(RPCA_RX_out,N,NN);

%% Post-detection algorithm: L2 norm of the column (RPCA_L2, the effect in HYDICE is better than RPCA_RX)
% E_spare = E_spare';
% RPCA_AD = zeros(1,N*NN);
% for i = 1: N*NN
%    RPCA_AD(:,i) = norm(E_spare(:,i));    
% end
% RPCA_AD = RPCA_AD';
% RPCA_RX_out = reshape(RPCA_AD,N,NN);
% 
% % ALM Perform RPCA calculations

end

