clear all
%% This test is trying level set method
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

%% level set scheme
tic
dy = 0.1;
y = 0:dy:2; Ny = length(y);
phi = zeros(Nx,Ny);
eta = zeros(Nx*Ny,Nt);
for j = 1:Ny
    phi(:,j) = y(j)*ones(Nx,1)-u0;
end
eta(:,1) = reshape(phi,[Nx*Ny,1]);

for i = 1:Nt-1
    Dx0 = (circshift(phi,[1,-1])-circshift(phi,[1,1]))/2/dx;
    Dy0 = (circshift(phi,[2,-1])-circshift(phi,[2,1]))/2/dy;
    Dyplus = (circshift(phi,[2,-1])-phi)/dy;
    Dyminus = (phi-circshift(phi,[2,1]))/dy;
    for j = 1:Ny
        phi(:,j) = phi(:,j)-dt/dx*y(j)*(phi(:,j)-circshift(phi(:,j),1));
    end
    eta(:,i+1) = reshape(phi(:,:),[Nx*Ny,1]);
end
toc

tic
%% Lagrangian level set
X = zeros(Nx,Ny);
Phi = zeros(Nx,Ny);
Xl = zeros(Nx*Ny,Nt);
Eta = zeros(Nx*Ny,Nt);
for j = 1:Ny
    Phi(:,j) = y(j)*ones(Nx,1)-u0;
end
for i = 1:Nt
    for j = 1:Ny
        X(:,j) = x+dt*(i-1)*y(j);
        Phi(:,j) = Phi(:,j);
    end
    Xl(:,i) = reshape(X,[Nx*Ny,1]);
    Eta(:,i) = reshape(Phi,[Nx*Ny,1]);
end
toc

%% now try DMD
tic
X = [Xl;Eta];
M = 250;
% [Z_k,Lambda_k,~,~,k] = DDMD_RRR(X(:,1:M),10^(-12));
%% DMD
X1 = X(:,1:M-1);
X2 = X(:,2:M);
[U,Sigma,V] = svd(X1,'econ');
index = find(diag(Sigma)<= sum(diag(Sigma))*1e-8);
k = min(index);
%k = 40;
U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*X2*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;%X2*V_k/Sigma_k*W;
Lambda_k = diag(D);
% 
%% DMD Spectra
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);
toc
%% Compute DMD Solution
x1 = X(:,1);
b = Z_k\x1;
time_dynamics = zeros(k,length(t));
for iter = 1:Nt
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Z_k*time_dynamics;
X_dmd = real(X_dmd);

space_dmd = zeros(Nx,Ny,Nt);
phi_dmd = zeros(Nx,Ny,Nt);
for i = 1:Nt
    space_dmd(:,:,i) = reshape(X_dmd(1:Nx*Ny,i),[Nx,Ny]);
    phi_dmd(:,:,i) = reshape(X_dmd(Nx*Ny+1:end,i),[Nx,Ny]);
end

Y = zeros(Nx,Ny);
for i = 1:Nx
    Y(i,:) = y;
end
[cc1] = contour(space_dmd(:,:,1)',Y',phi_dmd(:,:,1)',[0 0],'LineWidth',2);
[cc2] = contour(space_dmd(:,:,M)',Y',phi_dmd(:,:,M)',[0 0],'LineWidth',2);
[cc3] = contour(space_dmd(:,:,2*M)',Y',phi_dmd(:,:,2*M)',[0 0],'LineWidth',2);
[cc4] = contour(space_dmd(:,:,end)',Y',phi_dmd(:,:,end)',[0 0],'LineWidth',2);

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
hold on
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
plot(x,u(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(cc1(1,2:end),cc1(2,2:end),'k-.','LineWidth',1.2);
plot(cc2(1,2:end),cc2(2,2:end),'b-.','LineWidth',1.2);
plot(cc3(1,2:end),cc3(2,2:end),'g-.','LineWidth',1.2);
plot(cc4(1,2:end),cc4(2,2:end),'r-.','LineWidth',1.2);
axis([0 7 -0.1 2.1])
title('Level-Set DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
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

print -depsc2 level_set_dmd_burgers.eps


g_error = sqrt(sum((X-X_dmd).^2,1));

tau = 0*g_error;
for i = 2:M
    tau(i) = norm(Z_k\X(:,i)-Lambda_k*(Z_k\X(:,i-1)),2);
end

bd = g_error;
for i = M:Nt
    bd(i) = norm(Z_k,'fro')*(g_error(M)+(i-M)*max(tau(1:i)));
end

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
% xlabel({'t'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% ylabel({'error of g'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% 
% print -depsc2 error_bd_level_set_burgers.eps

figure
contour(space_dmd(:,:,1)',Y',phi_dmd(:,:,1)','-.','LineWidth',1.2);
hold on;
contour(space_dmd(:,:,1)',Y',phi_dmd(:,:,1)',[0 0],'r-','LineWidth',2);
axis([0 7 -0.1 2.1])
title('Level-Set Contours t =0','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'contours','0-level'},'FontSize',9,'Location','Southwest');

set(gca,'Units','normalized','Position',[.1 .1 .4 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'x'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'y'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 level_set1.eps

figure
contour(space_dmd(:,:,250)',Y',phi_dmd(:,:,250)','-.','LineWidth',1.2);
hold on;
contour(space_dmd(:,:,250)',Y',phi_dmd(:,:,250)',[0 0],'r-','LineWidth',2);
axis([0 7 -0.1 2.1])
title('Level-Set Contours t =0.25','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'contours','0-level'},'FontSize',9,'Location','Southwest');

set(gca,'Units','normalized','Position',[.1 .1 .4 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'x'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'y'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 level_set2.eps

figure
contour(space_dmd(:,:,500)',Y',phi_dmd(:,:,500)','-.','LineWidth',1.2);
hold on;
contour(space_dmd(:,:,500)',Y',phi_dmd(:,:,500)',[0 0],'r-','LineWidth',2);
axis([0 7 -0.1 2.1])
title('Level-Set Contours t =0.5','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'contours','0-level'},'FontSize',9,'Location','Southwest');

set(gca,'Units','normalized','Position',[.1 .1 .4 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'x'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'y'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 level_set3.eps

figure
contour(space_dmd(:,:,end)',Y',phi_dmd(:,:,end)','-.','LineWidth',1.2);
hold on;
contour(space_dmd(:,:,end)',Y',phi_dmd(:,:,end)',[0 0],'r-','LineWidth',2);
axis([0 7 -0.1 2.1])
title('Level-Set Contours t =1','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'contours','0-level'},'FontSize',9,'Location','Southwest');

set(gca,'Units','normalized','Position',[.1 .1 .4 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'x'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'y'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);

print -depsc2 level_set4.eps





% % subplot(1,2,2)
% plot(x,u(:,1),'k','LineWidth',2);
% hold on
% plot(x,u(:,M),'r','LineWidth',2);
% plot(x,u(:,2*M),'g','LineWidth',2);
% plot(x,u(:,end),'b','LineWidth',2);
% hold on;
% plot(x,u_dmd(:,1),'ko','LineWidth',2);
% plot(x,u_dmd(:,M),'ro','LineWidth',2);
% plot(x,u_dmd(:,2*M),'go','LineWidth',2);
% plot(x,u_dmd(:,end),'bo','LineWidth',2);
% xlabel('x');
% ylabel('u');
% title('reference solution vs Lagrangian-DMD solution');
% legend('reference t = 0','reference t = 0.3', 'reference t = 0.6','reference t = 1','DMD t = 0','DMD t = 0.3','DMD t = 0.6', 'DMD t = 1','Location','Best');
% 
% 
% 
% % % 
% % % subplot(1,2,2)
% % % error = sqrt(sum((X(Nx+1:end,:)-X_dmd(Nx+1:end,:)).^2,1));
% % % plot(t,error,'LineWidth',2);
% % % set(gca, 'YScale', 'log');
% % % xlabel('t');
% % % ylabel('error (L2 norm)');
% % % title('Lagrangian-DMD accuracy');
% % % 
