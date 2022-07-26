clear all
dt = 0.0001; t = 0:dt:1; Nt = length(t);
dx = 0.001; x = 0:dx:3; x = x'; Nx = length(x);


u = zeros(length(x),length(t));
u0 = ones(Nx,1);
u0(1:(Nx-1)/2) = 0.5;
u(:,1) = u0;

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

for i = 1:length(t)-1
    a_plus = max(u(:,i),0);
    a_minus = min(u(:,i),0);
    u_x_minus = (u(:,i)-circshift(u(:,i),1))/dx;
    u_x_plus = (circshift(u(:,i),-1)-u(:,i))/dx;
    u(:,i+1) = u(:,i)-dt*(a_plus.*u_x_minus+a_minus.*u_x_plus);
    u(1,i+1) = u(2,i+1);
    u(end,i+1) = u(end-1,i+1);
end

%% DMD
Y1 = zeros(Nx,Nt);
Y2 = zeros(Nx,Nt);
Y1(:,1) = x;
Y2(:,1) = u(:,1);
for i = 2:Nt
    Y1(:,i) = Y1(:,i-1)+dt*Y2(:,i-1);
    Y2(:,i) = Y2(:,i-1);
end
X = [Y1;Y2];

% M = 3000;
% plot(x,u(:,1),'k','LineWidth',2);
% hold on
% plot(x,u(:,M),'r','LineWidth',2);
% plot(x,u(:,2*M),'g','LineWidth',2);
% plot(x,u(:,end),'b','LineWidth',2);
% plot(Y1(:,1),Y2(:,1),'ko');
% plot(Y1(:,M),Y2(:,M),'ro');
% plot(Y1(:,2*M),Y2(:,2*M),'go');
% plot(Y1(:,end),Y2(:,end),'bo');


M = 3000;
[Z_k,Lambda_k,~,~,k] = DDMD_RRR(X(:,1:M),10^(-12));
% %%DMD
% X1 = X(:,1:M-1);
% X2 = X(:,2:M);
% [U,Sigma,V] = svd(X1,'econ');
% k = 20;
% U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
% Atilde = U_k'*X2*V_k/Sigma_k;
% [W,D] = eig(Atilde);
% Z_k = X2*V_k/Sigma_k*W;
% Lambda_k = diag(D);

%%DMD Spectra
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);

%% Compute DMD Solution
x1 = X(:,1);
b = Z_k\x1;
time_dynamics = zeros(k,length(t));
for iter = 1:Nt
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Z_k*time_dynamics;
X_dmd = real(X_dmd);

subplot(1,2,1)
plot(x,u(:,1),'k','LineWidth',2);
hold on
plot(x,u(:,M),'r','LineWidth',2);
plot(x,u(:,2*M),'g','LineWidth',2);
plot(x,u(:,end),'b','LineWidth',2);
hold on;
plot(X_dmd(1:100:Nx,1),X_dmd(Nx+1:100:end,1),'ko','LineWidth',2);
plot(X_dmd(1:100:Nx,M),X_dmd(Nx+1:100:end,M),'ro','LineWidth',2)
plot(X_dmd(1:100:Nx,2*M),X_dmd(Nx+1:100:end,2*M),'go','LineWidth',2);
plot(X_dmd(1:100:Nx,end),X_dmd(Nx+1:100:end,end),'bo','LineWidth',2);
xlabel('x');
ylabel('u');
title('reference solution vs Lagrangian-DMD solution');
legend('reference t = 0','reference t = 0.3', 'reference t = 0.6','reference t = 1','DMD t = 0','DMD t = 0.3','DMD t = 0.6', 'DMD t = 1','Location','Best');

subplot(1,2,2)
error = sqrt(sum((X(Nx+1:end,:)-X_dmd(Nx+1:end,:)).^2,1));
plot(t,error,'LineWidth',2);
set(gca, 'YScale', 'log');
xlabel('t');
ylabel('error (L2 norm)');
title('Lagrangian-DMD accuracy');

