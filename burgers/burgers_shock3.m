clear all
%% This test is trying relaxation scheme
dt = 0.0002; t = 0:dt:2; Nt = length(t);
dx = 0.002; x = 0:dx:3; x = x'; Nx = length(x);


% u = zeros(length(x),length(t));
% u0 = 0.8+0.5*exp(-(x-0.3).^2/0.01);
% u(:,1) = u0;

dt = 0.0001; t = 0:dt:1; Nt = length(t);
dx = 0.001; x = 0:dx:2; x = x'; Nx = length(x);


u = zeros(length(x),length(t));
u0 = 0.9*ones(Nx,1);
u0((Nx+1)/2:Nx) = 0.6;
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

%% relaxation scheme
a = 1;
ur = 0*u;
ur(:,1) = u0;
phi = 0*u;
eta = 0*u;
v = 0*u;
v0 = 0.5*u0.^2;
v(:,1) = v0;
phi (:,1) = u0+1/sqrt(a)*v0;
eta(:,1) = u0-1/sqrt(a)*v0;
eps = 1e-8;

for i = 1:Nt-1
    phi(:,i+1) = phi(:,i)-sqrt(a)*dt/dx*(phi(:,i)-circshift(phi(:,i),1));
    eta(:,i+1) = eta(:,i)+sqrt(a)*dt/dx*(circshift(eta(:,i),-1)-eta(:,i));
    phi(1,i+1) = phi(2,i+1);phi(end,i+1) = phi(end-1,i+1);
    eta(1,i+1) = eta(2,i+1);eta(end,i+1) = eta(end-1,i+1);
    ur(:,i+1) = 0.5*(phi(:,i+1)+eta(:,i+1));
    v(:,i+1) = sqrt(a)*0.5*(phi(:,i+1)-eta(:,i+1));
    v(:,i+1) = eps/(eps-dt)*(v(:,i+1)-dt/eps*0.5*ur(:,i+1).^2);
    phi(:,i+1) = u(:,i+1)+1/sqrt(a)*v(:,i+1);
    eta(:,i+1) = u(:,i+1)-1/sqrt(a)*v(:,i+1);
end

% %% check two solns consistent
% M = 3000;
% plot(x,u(:,1),'k','LineWidth',2);
% hold on
% plot(x,u(:,M),'r','LineWidth',2);
% plot(x,u(:,2*M),'g','LineWidth',2);
% plot(x,u(:,end),'b','LineWidth',2);
% plot(x,ur(:,1),'ko');
% plot(x,ur(:,M),'ro');
% plot(x,ur(:,2*M),'go');
% plot(x,ur(:,end),'bo');

%% now try DMD
X = [ur;v];
  M = 3000;
 [Z_k,Lambda_k,~,~,k] = DDMD_RRR(X(:,1:M),10^(-12));
% % %% DMD
% % X1 = X(:,1:M-1);
% % X2 = X(:,2:M);
% % [U,Sigma,V] = svd(X1,'econ');
% % k = rank(Sigma);
% % %k = 40;
% % U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
% % Atilde = U_k'*X2*V_k/Sigma_k;
% % [W,D] = eig(Atilde);
% % Z_k = X2*V_k/Sigma_k*W;
% % Lambda_k = diag(D);
% 
%% DMD Spectra
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

u_dmd = 0.5*(X_dmd(1:Nx,:)+X_dmd(Nx+1:2*Nx,:));
v_dmd = sqrt(a)*0.5*(X_dmd(1:Nx,:)-X_dmd(Nx+1:2*Nx,:));

 
% subplot(1,2,2)
plot(x,u(:,1),'k','LineWidth',2);
hold on
plot(x,u(:,M),'r','LineWidth',2);
plot(x,u(:,2*M),'g','LineWidth',2);
plot(x,u(:,end),'b','LineWidth',2);
hold on;
plot(x,u_dmd(:,1),'ko','LineWidth',2);
plot(x,u_dmd(:,M),'ro','LineWidth',2);
plot(x,u_dmd(:,2*M),'go','LineWidth',2);
plot(x,u_dmd(:,end),'bo','LineWidth',2);
xlabel('x');
ylabel('u');
title('reference solution vs Lagrangian-DMD solution');
legend('reference t = 0','reference t = 0.3', 'reference t = 0.6','reference t = 1','DMD t = 0','DMD t = 0.3','DMD t = 0.6', 'DMD t = 1','Location','Best');



% % 
% % subplot(1,2,2)
% % error = sqrt(sum((X(Nx+1:end,:)-X_dmd(Nx+1:end,:)).^2,1));
% % plot(t,error,'LineWidth',2);
% % set(gca, 'YScale', 'log');
% % xlabel('t');
% % ylabel('error (L2 norm)');
% % title('Lagrangian-DMD accuracy');
% % 
