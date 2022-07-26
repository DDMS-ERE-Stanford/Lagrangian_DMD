function [sigma,w] = svd_min(A)
[~,S,W] = svd(A);
sigma = S(end,end);
w = W(:,end);
end