clear all
tic
dt = 0.001; t = 0:dt:1; Nt = length(t);
dx = pi/1000; x = 0:dx:2*pi; x = x'; Nx = length(x);

u = zeros(length(x),length(t));
u0 = 1+sin(x);
u(:,1) = u0;

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

for i = 1:length(t)-1
    F = (0.5*u(:,i).^2+0.5*(circshift(u(:,i),-1)).^2)/2-...
        dx/2/dt*(circshift(u(:,i),-1)-u(:,i)); %F_{i+1/2}
    u(:,i+1) = u(:,i)-dt/dx*(F-circshift(F,1));
end
toc
%% DMD
tic
Y1 = zeros(Nx,Nt);
Y2 = zeros(Nx,Nt);
Y1(:,1) = x;
Y2(:,1) = u(:,1);
for i = 2:Nt
%     x_hat = Y1(:,i-1)+dt/2*Y2(:,i-1);
%     Y1(:,i) = Y1(:,i-1)+dt*interp1(x,0.5*(u(:,i)+u(:,i-1)),x_hat,'spline','extrap');
    Y2(:,i) = Y2(:,i-1);
    Y1(:,i) = Y1(:,i-1)+dt*Y2(:,i);
end
X = [Y1;Y2];
toc
M = 250;

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
Z_k = U_k*W;%X2*V_k/Sigma_k*W;
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

%%POD
tic
[U,Sigma,V] = svd(X,'econ');
U_k = U(:,1:k);
y_hat = zeros(k,Nt);
y_pod = zeros(2*Nx,Nt);
y_pod(:,1) =[x;u0];
y_hat(:,1) = U_k\[x;u0];
% A2_hat = U_k'*[diag(ones(Nx,1)),zeros(Nx,Nx);zeros(Nx,Nx),A2]*U_k;

for i = 1:length(t)-1
    y_hat(:,i+1) = y_hat(:,i)+dt*U_k'*[y_pod(Nx+1:end,i);zeros(Nx,1)];
    y_pod(:,i+1) = U_k*y_hat(:,i+1);
end
    u_pod = y_pod(Nx+1:end,:);
    x_pod = y_pod(1:Nx,:);
toc

dt = 0.001; t = 0:dt:1; Nt = length(t);
dx = pi/1000; x = 0:dx:4*pi; x = x';

u = zeros(length(x),length(t));
u0 = 1+sin(x);
u(:,1) = u0;

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

for i = 1:length(t)-1
    F = (0.5*u(:,i).^2+0.5*(circshift(u(:,i),-1)).^2)/2-...
        dx/2/dt*(circshift(u(:,i),-1)-u(:,i)); %F_{i+1/2}
    u(:,i+1) = u(:,i)-dt/dx*(F-circshift(F,1));
end

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


axis([0 7 -0.1 2.1])
title('DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'ref $t = 0$','ref $t = 0.25$',...
    'ref $t = 0.5$','ref $t = 1$',...
    'DMD $t = 0$','DMD $t = 0.25$','DMD $t = 0.5$','DMD $t = 1$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Bestoutside');
legend('boxoff');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 dmd_burgers.eps

figure
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(x_pod(:,1),u_pod(:,1),'k-.','LineWidth',1.2);
plot(x_pod(:,250),u_pod(:,250),'b-.','LineWidth',1.2);
plot(x_pod(:,500),u_pod(:,500),'g-.','LineWidth',1.2);
plot(x_pod(:,end),u_pod(:,end),'r-.','LineWidth',1.2);

axis([0 7 -0.1 2.1])
title('POD vs. Reference Solution','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'ref $t = 0$','ref $t = 0.25$',...
    'ref $t = 0.5$','ref $t = 1$',...
    'POD $t = 0$','POD $t = 0.25$','POD $t = 0.5$','POD $t = 1$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Bestoutside');
legend('boxoff');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 pod_burgers.eps


% ref_soln = 0*u_pod;
% for i =1:Nt
%     ref_soln(:,i) = interp1(x,u(:,i),X_dmd(1:Nx,i),'linear');
% end
% 
% error_dmd = sqrt(sum((ref_soln-X_dmd(Nx+1:end,:)).^2,1));
% for i =1:Nt
%     ref_soln(:,i) = interp1(x,u(:,i),x_pod(1:Nx,i),'linear');
% end
% error_pod = sqrt(sum((ref_soln-u_pod).^2,1));
error_dmd = sqrt(sum((Y2-X_dmd(Nx+1:end,:)).^2,1));
error_pod = sqrt(sum((Y2-u_pod).^2,1));

figure
plot(t,error_dmd,'LineWidth',2);
hold on;
plot(t,error_pod,'-.','LineWidth',2);

title('DMD and POD accuracy','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'DMD','POD'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');
legend('boxoff');

set(gca,'YScale', 'log','Position',[.1 .1 0.4 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$t$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'error of $\bf u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 error_burgers.eps

g_error = sqrt(sum(([Y1;Y2]-X_dmd).^2,1));
tau = 0*g_error;
for i = 2:M
    tau(i) = norm(Z_k\[Y1(:,i);Y2(:,i)]-Lambda_k*(Z_k\[Y1(:,i-1);Y2(:,i-1)]),2);
end

bd = g_error;
for i = M:Nt
    bd(i) = norm(Z_k,'fro')*(g_error(M)+(i-M)*max(tau));
end

figure
plot(t(M:end),g_error(M:end),'k','LineWidth',2);
hold on;
plot(t(M:end),bd(M:end),'k-.','LineWidth',2);
set(gca,'YScale', 'log','Position',[.1 .1 .4 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
title('DMD error estimate','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'DMD observable error','error bound'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');
legend('boxoff');
xlabel({'$t$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'error of $\bf g(\bf u)$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 error_bd_burgers.eps
