function [k,var,U,E,l] = wic(snapshots)

%% Weighted Information Criteria
% Input:
% Snapshots: Snapshots of data. The covariance is calculated inside
% k = Optimal model order
% var = Maximum Likelihood Estimation of variance
% U = Left singular vectors
% E = Singular Values
% l = Eigenvalues of covariance matrix of snapshots
[p n] = size(snapshots);
[U E V] = svd(snapshots);


l = diag(E).^2/(2*n);
wic = zeros(1,p-1);
for k=1:p-1
    M = k*(2*p-k);
    N = n;
    W = ((2*M*N/(N-M-1))^2 + (M*log(N))^2)/(2*N*M/(N-M-1) + M*log(N));
    %aic(k) = 2*k*(2*p-k) -n*(sum(log(l(k+1:p))) - (p-k)*log(mean(l(k+1:p))));
    %aic(k) = -2*n*sum(log(l(k+1:p))) + 2*n*(p-k)*log(mean(l(k+1:p))) + 2*k*(2*p-k);
    wic(k) = -2*n*sum(log(l(k+1:p))) + 2*n*(p-k)*log(mean(l(k+1:p))) + W;
end

%figure(2)
%plot(aic)
[mn,k] = min(wic);
var = sum(l(k+1:p))/(p-k);
end