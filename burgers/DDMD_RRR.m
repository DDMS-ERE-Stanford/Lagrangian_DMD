function [Z_k,Lambda_k,r_k,rho_k,k] = DDMD_RRR(X,eps)
D_x = diag(sqrt(sum(X(:,1:end-1).^2,1)));
X1 = X(:,1:end-1)*D_x;
Y1 = X(:,2:end)*D_x;
[U,Sigma,V] = svd(X1);
index = find(diag(Sigma)>= Sigma(1,1)*eps);
k = max(index);
U_k = U(:,1:k);
V_k = V(:,1:k);
Sigma_k = Sigma(1:k,1:k);
B_k = Y1*(V_k/Sigma_k);
[~,R] = qr([U_k,B_k]);
S_k = diag(diag(R(1:k,1:k)))*R(1:k,k+1:2*k);
Lambda_k = eig(S_k);
W_k = zeros(k,k);
r_k = zeros(k,1);
rho_k = zeros(k,1);
for i = 1:k
    [sigma,w] = svd_min([R(1:k,k+1:2*k)-Lambda_k(i)*R(1:k,1:k);R(k+1:2*k,k+1:2*k)]);
    W_k(:,i) = w;
    r_k(i) = sigma;
    rho_k(i) = w'*S_k*w;
end
Z_k = U_k*W_k;
end