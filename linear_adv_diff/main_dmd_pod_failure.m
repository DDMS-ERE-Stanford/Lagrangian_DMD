clear all
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
A2 = Id-dt*diffusion*D2;

for i = 1:length(t)-1
    u(:,i+1) = A2\(A1*u(:,i));
end

% DMD
M = 100;
X = u(:,1:M);
% X = u(:,1:5:end);

[U,Sigma,V] = svd(X,'econ');
figure
width = 3;     % Width in inches
height = 3;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

plot(diag(Sigma),'k-','LineWidth',1);
set(gca, 'YScale', 'log')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',8)
xlabel({'index'},'FontUnits','points','interpreter','latex',...
    'FontSize',8);
ylabel({'singular value'},'FontUnits','points','interpreter','latex',...
    'FontSize',8);
% print('svd2.eps','-depsc2','-r300');



% 
% 
% 
% 
% 
% 
% 
% 
% 
% X1 = X(:,1:M-1);
% X2 = X(:,2:M);
% [U,Sigma,V] = svd(X1,'econ');
% % index = find(diag(Sigma)<= sum(diag(Sigma))*1e-14);
% % k = min(index);
% k = 20;
% U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
% Atilde = U_k'*X2*V_k/Sigma_k;
% [W,D] = eig(Atilde);
% Z_k = U_k*W;
% Lambda_k = diag(D);
% 
% %DMD Spectra
% omega = log(Lambda_k)/dt;
% Lambda_k = diag(Lambda_k);
% 
% % Compute DMD Solution
% x1 = X(:,1);
% b = Z_k\x1;
% time_dynamics = zeros(k,length(t));
% for iter = 1:Nt
%     time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
% end
% X_dmd = Z_k*time_dynamics;
% X_dmd = real(X_dmd);
% 
% %POD
% [U,Sigma,V] = svd(X,'econ');
% U_k = U(:,1:k);
% u_hat = zeros(k,Nt);
% u_hat(:,1) = U_k\u0;
% u_pod = 0*u;
% u_pod(:,1) = u0;
% A1_hat = U_k'*A1*U_k;
% A2_hat = U_k'*A2*U_k;
% 
% for i = 1:length(t)-1
%     u_hat(:,i+1) = A2_hat\(A1_hat*u_hat(:,i));
%     u_pod(:,i+1) = U_k*u_hat(:,i+1);
% end
% 
% 
% figure
% plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
% hold on
% plot(x,u(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
% plot(x,u(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
% plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
% plot(x,X_dmd(:,1),'k-.','LineWidth',1.2);
% plot(x,X_dmd(:,250),'b-.','LineWidth',1.2);
% plot(x,X_dmd(:,500),'g-.','LineWidth',1.2);
% plot(x,X_dmd(:,end),'r-.','LineWidth',1.2);
% 
% 
% % axis([x(1) x(end) -0.3 1])
% title('DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
% legend({'ref $t = 0$','ref $t = 0.25$',...
%     'ref $t = 0.5$','ref $t = 1$',...
%     'DMD $t = 0$','DMD $t = 0.25$','DMD $t = 0.5$','DMD $t = 1$'},...
%     'FontUnits','points','interpreter','latex',...
%     'FontSize',9,'Location','northeast');
% 
% set(gca,'Units','normalized','Position',[.1 .1 .4 .7],...
%     'FontUnits','points','FontWeight','normal','FontSize',9)
% xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% legend('boxoff');
% print -depsc2 dmd2.eps
% 
% figure
% plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
% hold on
% plot(x,u(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
% plot(x,u(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
% plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
% plot(x,u_pod(:,1),'k-.','LineWidth',1.2);
% plot(x,u_pod(:,250),'b-.','LineWidth',1.2);
% plot(x,u_pod(:,500),'g-.','LineWidth',1.2);
% plot(x,u_pod(:,end),'r-.','LineWidth',1.2);
% 
% 
% % axis([x(1) x(end) -0.3 1])
% title('POD vs. Reference Solution','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
% legend({'ref $t = 0$','ref $t = 0.25$',...
%     'ref $t = 0.5$','ref $t = 1$',...
%     'POD $t = 0$','POD $t = 0.25$','POD $t = 0.5$','POD $t = 1$'},...
%     'FontUnits','points','interpreter','latex',...
%     'FontSize',9,'Location','northeast');
% 
% set(gca,'Units','normalized','Position',[.1 .1 .4 .7],...
%     'FontUnits','points','FontWeight','normal','FontSize',9)
% xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
%     'FontSize',9);
% legend('boxoff');
% print -depsc2 pod2.eps
% 
% % error_dmd = sqrt(sum((u-X_dmd).^2,1));
% % error_pod = sqrt(sum((u-u_pod).^2,1));
% % figure
% % plot(t,error_dmd,'LineWidth',2);
% % hold on;
% % plot(t,error_pod,'LineWidth',2);
% % 
% % title('DMD and POD accuracy','FontUnits','points','interpreter','latex',...
% %     'FontSize',10)
% % legend({'DMD','POD'},...
% %     'FontUnits','points','interpreter','latex',...
% %     'FontSize',9,'Location','Best');
% % 
% % set(gca,'YScale', 'log','Position',[.1 .1 .8 .4],...
% %     'FontUnits','points','FontWeight','normal','FontSize',9)
% % xlabel({'t'},'FontUnits','points','interpreter','latex',...
% %     'FontSize',9);
% % ylabel({'error'},'FontUnits','points','interpreter','latex',...
% %     'FontSize',9);
% % 
% % print -depsc2 error1.eps
% % 
% % figure
% % hold on;
% % plot(x,Z_k(:,1),'k-','LineWidth',1);
% % plot(x,Z_k(:,10),'b-.','LineWidth',1);
% % plot(x,Z_k(:,20),'g--','LineWidth',1);
% % 
% % 
% % axis([x(1) x(end) -0.1 0.1])
% % title('DMD modes','FontUnits','points','interpreter','latex',...
% %     'FontSize',10)
% % legend({'DMD mode1','DMD mode10','DMD mode20'},...
% %     'FontUnits','points','interpreter','latex',...
% %     'FontSize',9,'Location','northeast');
% % 
% % set(gca,'Units','normalized','Position',[.1 .1 .4 .7],...
% %     'FontUnits','points','FontWeight','normal','FontSize',9)
% % xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
% %     'FontSize',9);
% % ylabel({'$\phi$'},'FontUnits','points','interpreter','latex',...
% %     'FontSize',9);
% % legend('boxoff');
% % print -depsc2 dmd_modes.eps
% % 
% % figure
% % hold on;
% % plot(x,U_k(:,1),'k-','LineWidth',1);
% % plot(x,U_k(:,10),'b-.','LineWidth',1);
% % plot(x,U_k(:,20),'g--','LineWidth',1);
% % 
% % 
% % axis([x(1) x(end) -0.1 0.1])
% % title('POD modes','FontUnits','points','interpreter','latex',...
% %     'FontSize',10)
% % legend({'POD mode1','POD mode10','POD mode20'},...
% %     'FontUnits','points','interpreter','latex',...
% %     'FontSize',9,'Location','northeast');
% % 
% % set(gca,'Units','normalized','Position',[.1 .1 .4 .7],...
% %     'FontUnits','points','FontWeight','normal','FontSize',9)
% % xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
% %     'FontSize',9);
% % ylabel({'$\phi$'},'FontUnits','points','interpreter','latex',...
% %     'FontSize',9);
% % legend('boxoff');
% % print -depsc2 pod_modes.eps
% % 
% 
% 
