clear all;
tic
dt = 0.001; t = 0:dt:1; Nt = length(t);
dx = 0.001; x = 0:dx:2; x = x'; Nx = length(x);

u = zeros(length(x),length(t));
u0 = 0.5*exp(-(x-0.3).^2/0.05^2);
u(:,1) = u0;

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

advection = 1;
diffusion = 0;


vtmp = ones(Nx-1,1);
D1 = (diag(ones(Nx,1))-diag(vtmp,-1))/dx;
D1(1,:) = 0; D1(end,:)=0;
D2 = (diag(-2*ones(Nx,1))+diag(vtmp,1)+diag(vtmp,-1))/(dx^2);   % second derivative matrix
D2(1,:) = 0; D2(end,:) =0;
Id = diag(ones(Nx,1));     % identity matrix
A1 = Id-dt*advection*D1;

for i = 1:length(t)-1
    u(:,i+1) = A1*u(:,i);
end
toc


%% DMD
tic
M = 250;
Xl = zeros(Nx,Nt);
Ul = zeros(Nx,Nt);
for i = 1:Nt
    Xl(:,i) = x+dt*(i-1)*advection;
    Ul(:,i) = u(:,1);
end
X = [Xl(:,1:M);Ul(:,1:M)];
toc

%%DMD
tic
X1 = X(:,1:M-1);
X2 = X(:,2:M);
[U,Sigma,V] = svd(X1,'econ');
index = find(diag(Sigma)<= sum(diag(Sigma))*1e-8);
k = min(index);
% k = 5;
U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*X2*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;
Lambda_k = diag(D);
toc

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

figure
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(X_dmd(1:Nx,1),X_dmd(Nx+1:end,1),'k-.','LineWidth',1.2);
plot(X_dmd(1:Nx,250),X_dmd(Nx+1:end,250),'b-.','LineWidth',1.2);
plot(X_dmd(1:Nx,500),X_dmd(Nx+1:end,500),'g-.','LineWidth',1.2);
plot(X_dmd(1:Nx,end),X_dmd(Nx+1:end,end),'r-.','LineWidth',1.2);


axis([x(1) x(end) 0 1])
title('DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'ref $t = 0$','ref $t = 0.25$',...
    'ref $t = 0.5$','ref $t = 1$',...
    'DMD $t = 0$','DMD $t = 0.25$','DMD $t = 0.5$','DMD $t = 1$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','bestoutside');
legend('boxoff');

set(gca,'Units','normalized','Position',[.1 .1 .4 .7],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 dmd_advection.eps

% %%POD
% tic
% [U,Sigma,V] = svd(X,'econ');
% U_k = U(:,1:k);
% y_hat = zeros(k,Nt);
% y_pod = zeros(2*Nx,Nt);
% y_pod(:,1) =[x;u0];
% y_hat(:,1) = U_k\[x;u0];
% 
% for i = 1:length(t)-1
%     y_hat(:,i+1) = y_hat(:,i)+dt*U_k'*[advection*ones(Nx,1);zeros(Nx,1)];
%     y_pod(:,i+1) = U_k*y_hat(:,i+1);
% end
%     u_pod = y_pod(Nx+1:end,:);
%     x_pod = y_pod(1:Nx,:);
% toc
% 
% figure
% plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
% hold on
% plot(x,u(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
% plot(x,u(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
% plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
% plot(x_pod(:,1),u_pod(:,1),'k-.','LineWidth',1.2);
% plot(x_pod(:,250),u_pod(:,250),'b-.','LineWidth',1.2);
% plot(x_pod(:,500),u_pod(:,500),'g-.','LineWidth',1.2);
% plot(x_pod(:,end),u_pod(:,end),'r-.','LineWidth',1.2);
% 
% axis([x(1) x(end) 0 1])
% title('POD vs. Reference Solution','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
% legend({'ref $t = 0$','ref $t = 0.25$',...
%     'ref $t = 0.5$','ref $t = 1$',...
%     'POD $t = 0$','POD $t = 0.25$','POD $t = 0.5$','POD $t = 1$'},...
%     'FontUnits','points','interpreter','latex',...
%     'FontSize',9,'Location','northeast');
% legend('boxoff');
% 
% set(gca,'Units','normalized','Position',[.1 .1 .4 .7],...
%     'FontUnits','points','FontWeight','normal','FontSize',9)
% xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% 
% print -depsc2 pod_advection.eps
% 
% 
% error_dmd = sqrt(sum((u0-X_dmd(Nx+1:end,:)).^2,1));
% error_pod = sqrt(sum((u0-u_pod).^2,1));
% figure
% plot(t,error_dmd,'LineWidth',2);
% hold on;
% plot(t,error_pod,'-.','LineWidth',2);
% 
% title('DMD and POD accuracy','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
% legend({'DMD','POD'},...
%     'FontUnits','points','interpreter','latex',...
%     'FontSize',9,'Location','Best');
% legend('boxoff');
% 
% set(gca,'YScale', 'log','Position',[.1 .1 .4 .4],...
%     'FontUnits','points','FontWeight','normal','FontSize',9)
% xlabel({'$t$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% ylabel({'error of $\bf{u}$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% 
% print -depsc2 error_advection.eps
% 
% g_error = sqrt(sum(([Xl;Ul]-X_dmd).^2,1));
% tau = 0*g_error;
% for i = 2:Nt
%     tau(i) = norm(Z_k\[Xl(:,i);Ul(:,i)]-Lambda_k*(Z_k\[Xl(:,i-1);Ul(:,i-1)]),2);
% end
% 
% bd = g_error;
% for i = M:Nt
%     bd(i) = norm(Z_k,'fro')*(g_error(M)+(i-M)*max(tau));
% end
% 
% figure
% plot(t(M:end),g_error(M:end),'k','LineWidth',2);
% hold on;
% plot(t(M:end),bd(M:end),'k-.','LineWidth',2);
% set(gca,'YScale', 'log','Position',[.1 .1 .4 .4],...
%     'FontUnits','points','FontWeight','normal','FontSize',9)
% title('DMD error estimate','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
% legend({'DMD observable error','error bound'},...
%     'FontUnits','points','interpreter','latex',...
%     'FontSize',9,'Location','Best');
% legend('boxoff');
% xlabel({'$t$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% ylabel({'error of $\bf{g}(\bf{u})$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% 
% print -depsc2 error_bd_advection.eps
% 
% 
% 
% 
