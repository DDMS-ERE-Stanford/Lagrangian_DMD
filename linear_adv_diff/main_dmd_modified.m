clear all
dt = 0.01; t = 0:dt:2; Nt = length(t);
dx = 0.001; x = -2:dx:2; x = x'; Nx = length(x);

u = zeros(length(x),length(t));
advection = 0;
diffusion = 0.01;
% u0 = 0.5*exp(-(x-1).^2/0.05^2);

%% analytic solution u = 1/sqrt(4*pi*diffusion*(t+1))*exp(-x.^2/(4*diffusion*(t+1)));
u = zeros(length(x),length(t));
for i = 1:length(t)
    u(:,i) = 1/sqrt(4*pi*diffusion*(t(i)+1))*exp(-x.^2/(4*diffusion*(t(i)+1)));
end

% u1 = 0*u;
% u0 = 1/sqrt(4*pi*diffusion)*exp(-x.^2/4/diffusion);
% u1(:,1) = u0;
% % check CFL condition
% CFL = max(abs(u(:,1)))*dt/dx;
% fprintf('CFL number = %7.3f\n',CFL);
% 
% vtmp = ones(Nx-1,1);
% D1 = (diag(ones(Nx,1))-diag(vtmp,-1))/dx;
% D1(1,:) = 0; D1(end,:)=0;
% D2 = (diag(-2*ones(Nx,1))+diag(vtmp,1)+diag(vtmp,-1))/(dx^2);   % second derivative matrix
% D2(1,:) = 0; D2(end,:) =0;
% Id = diag(ones(Nx,1));     % identity matrix
% A1 = Id-dt*advection*D1;
% A2 = Id-dt*diffusion*D2;
% 
% for i = 1:length(t)-1
%     u1(:,i+1) = A2\(A1*u1(:,i));
% end

% %% DMD
% M = 250;
% X = u(:,1:M);
% X1 = X(:,1:M-1);
% X2 = X(:,2:M);
% [U,Sigma,V] = svd(X1,'econ');
% index = find(diag(Sigma)<= sum(diag(Sigma))*1e-4);
% k = min(index);
% % k = 10;
% U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
% Atilde = U_k'*X2*V_k/Sigma_k;
% [W,D] = eig(Atilde);
% % Z_k = U_k*W;%X2*V_k/Sigma_k*W;
% Z_k = X2*V_k/Sigma_k*W;
% Lambda_k = diag(D);
% 
% %% DMD Spectra
% omega = log(Lambda_k)/dt;
% Lambda_k = diag(Lambda_k);
% 
% %% Compute DMD Solution
% x1 = X(:,1);
% b = Z_k\x1;
% X_dmd = 0*u;
% X_dmd(:,1) = u0;
% for iter = 2:Nt
%     X_dmd(:,iter) = Z_k*Lambda_k.^(iter-1)*b;
% end
% X_dmd = real(X_dmd);

%% DMD
% M = 250;
M = Nt;
X = u(:,1:M);
X1 = X(:,1:M-1);
X2 = X(:,2:M);
[U,Sigma,V] = svd(X1,'econ');
index = find(diag(Sigma)<= sum(diag(Sigma))*1e-8);
k = min(index);
% k = 10;
U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*(X2-X1)/dt*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;%X2*V_k/Sigma_k*W;
% Z_k = X2*V_k/Sigma_k*W;
Lambda_k = dt*diag(D)+ones(k,1);

%% DMD Spectra
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);

%% Compute DMD Solution
x1 = X1(:,1);
b = Z_k\x1;
X_dmd = 0*u;
X_dmd(:,1) = u(:,1);
for iter = 2:Nt
    X_dmd(:,iter) = Z_k*Lambda_k.^(iter-1)*b;
end
X_dmd = real(X_dmd);

% figure
% plot(x,u(:,1),'LineWidth',6,'Color',[0.75,0.75,0.75]);
% hold on
% plot(x,u(:,250),'LineWidth',6,'Color',[0.6,0.8,1]);
% % plot(x,u(:,500),'LineWidth',6,'Color',[0.8,1,0.8]);
% % plot(x,u(:,end),'LineWidth',6,'Color',[1,0.6,0.6]);
% plot(x,X_dmd(:,1),'k-.','LineWidth',1.2);
% plot(x,X_dmd(:,250),'b-.','LineWidth',1.2);
% % plot(x,X_dmd(:,500),'g-.','LineWidth',1.2);
% % plot(x,X_dmd(:,end),'r-.','LineWidth',1.2);
% 
% xlabel('x');
% ylabel('u');
% title('DMD vs. Reference Solution')
% legend('reference soln t = 0','reference soln t = 0.25','reference soln t = 0.5','reference soln t = 1','DMD soln t = 0','DMD soln t = 0.25','DMD soln t = 0.5','DMD soln t = 1','Location','Best')

figure
hold on
error_dmd = sqrt(sum((u-X_dmd).^2,1));
plot(t,error_dmd,'LineWidth',2);
set(gca, 'YScale', 'log');
xlabel('t');
ylabel('error (L2 norm)');
title('DMD and POD accuracy');
legend('DMD','Location','Best')

