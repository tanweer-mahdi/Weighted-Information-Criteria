function [k,var,U,E,l] = aic(R,n)

%% Minimum Description Length (MDL) based Moder Order Selection (also known as Bayesian Information Criterion)
% Input:
% R = Covariance matrix of data snapshots
% n = Number of data snapshots
% Output:
% k = Optimal Model Order
% var = Maximum Likelihood Estimation of variance
[p ~] = size(R);
[eigvec eigvals] = eig(R);

l = sort(diag(eigvals),'descend');
aic = zeros(1,p-1);
for k=1:p-1
    %aic(k) = 2*k*(2*p-k) -n*(sum(log(l(k+1:p))) - (p-k)*log(mean(l(k+1:p))));
    %aic(k) = -2*n*sum(log(l(k+1:p))) + 2*n*(p-k)*log(mean(l(k+1:p))) + 2*k*(2*p-k);
    aic(k) = -n*sum(log(l(k+1:p))) + n*(p-k)*log(mean(l(k+1:p))) + 0.5*k*(2*p-k)*log(n);
end

%figure(2)
%plot(aic)
[mn,k] = min(aic);
var = sum(l(k+1:p))/(p-k);
end