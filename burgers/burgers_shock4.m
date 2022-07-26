%clear all
%% This test is trying level set method
dx = 0.02; x = 0:dx:2; x = x'; Nx = length(x);
dt = dx^3; t = 0:dt:0.5; Nt = length(t);

u = zeros(length(x),length(t));
u0 = 0.8+0.5*exp(-(x-0.3).^2/0.01);
u(:,1) = u0;

% u = zeros(length(x),length(t));
% u0 = ones(Nx,1);
% u0((Nx+1)/2:Nx) = 0.1;
% u(:,1) = u0;


% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

for i = 1:Nt-1
    a_plus = max(u(:,i),0);
    a_minus = min(u(:,i),0);
    u_x_minus = (u(:,i)-circshift(u(:,i),1))/dx;
    u_x_plus = (circshift(u(:,i),-1)-u(:,i))/dx;
    u(:,i+1) = u(:,i)-dt*(a_plus.*u_x_minus+a_minus.*u_x_plus);
    u(1,i+1) = u(2,i+1);
    u(end,i+1) = u(end-1,i+1);
end


%% level set scheme
dy = 0.1;
y = 0.7:dy:1.5; Ny = length(y);
phi = zeros(Nx,Ny);
eta = zeros(Nx*Ny,Nt);
for j = 1:Ny
    phi(:,j) = y(j)*ones(Nx,1)-u0;
end
eta(:,1) = reshape(phi,[Nx*Ny,1]);

m = 0.2;

for i = 1:Nt-1
    Dx0 = (circshift(phi,[-1,0])-circshift(phi,[1,0]))/2/dx;
    Dy0 = (circshift(phi,[0,-1])-circshift(phi,[0,1]))/2/dy;
    Dx0(1,:) = 0; Dx0(end,:) = 0; Dx0(:,1) = 0; Dx0(:,end) = 0;
    Dy0(:,1) = 0; Dy0(:,end) = 0; Dy0(1,:) = 0; Dy0(end,:) = 0;
    Dyplus = (circshift(phi,[0,-1])-phi)/dy;
    Dyplus(:,end) = 0;
    Dyminus = (phi-circshift(phi,[0,1]))/dy;
    Dyminus(:,1) = 0;
    for j = 1:Ny
        phi(:,j) = phi(:,j)-dt/dx*y(j)*(phi(:,j)-[phi(1,j);phi(1:end-1,j)])...
                +dt*m*sqrt(Dx0(:,j).^2+Dy0(:,j).^2).*(tanh(1/dy*Dyplus(:,j))-tanh(1/dy*Dyminus(:,j)))/dy;
    end
    eta(:,i+1) = reshape(phi(:,:),[Nx*Ny,1]);
end

        
%% check two solns consistent
figure
hold on
plot(x,u(:,end),'o','LineWidth',2);
contour(x,y,phi',[0 0],'LineWidth',2);

cc = contour(x,y,phi',[0 0]); % Overlay contour line


% %% Lagrangian level set
% X = zeros(Nx,Ny);
% Phi = zeros(Nx,Ny);
% Xl = zeros(Nx*Ny,Nt);
% Eta = zeros(Nx*Ny,Nt);
% for j = 1:Ny
%     Phi(:,j) = y(j)*ones(Nx,1)-u0;
% end
% for i = 1:Nt
%     Dx0 = (circshift(Phi,[-1,0])-circshift(Phi,[1,0]))/2/dx;
%     Dy0 = (circshift(Phi,[0,-1])-circshift(Phi,[0,1]))/2/dy;
%     Dx0(1,:) = 0; Dx0(end,:) = 0; Dx0(:,1) = 0; Dx0(:,end) = 0;
%     Dy0(:,1) = 0; Dy0(:,end) = 0; Dy0(1,:) = 0; Dy0(end,:) = 0;
%     Dyplus = (circshift(Phi,[0,-1])-Phi)/dy;
%     Dyplus(:,end) = 0;
%     Dyminus = (Phi-circshift(Phi,[0,1]))/dy;
%     Dyminus(:,1) = 0;
%     Dxplus = (circshift(Phi,[-1,0])-Phi)/dx;
%     Dxplus(end,:) = 0;
%     Dxminus = (Phi-circshift(Phi,[1,0]))/dx;
%     Dxminus(1,:) = 0;
%     for j = 1:Ny
%         X(:,j) = x+dt*(i-1)*y(j);
%         Phi(:,j) = Phi(:,j)...
%             +dt*m*sqrt(Dx0(:,j).^2+(Dx0(:,j)*t(i)+Dy0(:,j)).^2).*...
%             ((tanh(1/dy*(Dxplus(:,j)*t(i)+Dyplus(:,j)))-tanh(1/dy*(Dxminus(:,j)*t(i)+Dyminus(:,j))))/dy...
%             +t(i)*(tanh(1/dy*(Dxplus(:,j)*t(i)+Dyplus(:,j)))-tanh(1/dy*(Dxminus(:,j)*t(i)+Dyminus(:,j))))/dx);
%     end
%     Xl(:,i) = reshape(X,[Nx*Ny,1]);
%     Eta(:,i) = reshape(Phi,[Nx*Ny,1]);
% end
% 
% hold on
% plot(x,u(:,end),'o','LineWidth',2);
% Y = 0*X;
% for i = 1:Nx
%     Y(i,:) = y;
% end
% contour(X',Y',Phi',[0 0],'LineWidth',2);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %% now try DMD
% % X = [ur;v];
% %   M = 3000;
% %  [Z_k,Lambda_k,~,~,k] = DDMD_RRR(X(:,1:M),10^(-12));
% % % % %% DMD
% % % % X1 = X(:,1:M-1);
% % % % X2 = X(:,2:M);
% % % % [U,Sigma,V] = svd(X1,'econ');
% % % % k = rank(Sigma);
% % % % %k = 40;
% % % % U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
% % % % Atilde = U_k'*X2*V_k/Sigma_k;
% % % % [W,D] = eig(Atilde);
% % % % Z_k = X2*V_k/Sigma_k*W;
% % % % Lambda_k = diag(D);
% % % 
% % %% DMD Spectra
% % omega = log(Lambda_k)/dt;
% % Lambda_k = diag(Lambda_k);
% % 
% % %% Compute DMD Solution
% % x1 = X(:,1);
% % b = Z_k\x1;
% % time_dynamics = zeros(k,length(t));
% % for iter = 1:Nt
% %     time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
% % end
% % X_dmd = Z_k*time_dynamics;
% % X_dmd = real(X_dmd);
% % 
% % u_dmd = 0.5*(X_dmd(1:Nx,:)+X_dmd(Nx+1:2*Nx,:));
% % v_dmd = sqrt(a)*0.5*(X_dmd(1:Nx,:)-X_dmd(Nx+1:2*Nx,:));
% 
