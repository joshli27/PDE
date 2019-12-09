%% 1D Burgers' equation 
clear all
close all

%% set parameters
N          = 256;
deltax     = 1/N;
xmesh      = 0:deltax:1;
couxx      = 0;
ax         = couxx*ones(1,N+1);     % xmesh;
f          = @(x) sin(4*pi*x);      % initial condition
df         = @(x) 4*pi*cos(4*pi*x); % derivative of initial
Bleft      = @(x) 0; 
Bright     = @(x) 0; 
u0         = f(xmesh);

deltat     = 0.001; %deltax^2/2;
T          = 0.05;
tmesh      = 0:deltat:T;
LenT       = length(tmesh);

IfNormCol  = 1;

noise      = 10;
lambda     = 500;

u          = zeros(LenT,N+1);
u(1,:)     = u0;


% %% exact solution of pde
% Newton       = struct();
% Newton.Error = zeros(LenT,N+1);
% Newton.Imax  = zeros(LenT,N+1);
% Fzero        = struct();
% Fzero.Error  = zeros(LenT,N+1);
% for itert = 2 : LenT
%     for iterx = 2 : N
%         % Newton method
%         [ x0 , ex0 ]              = newton(@(x) f(x)*tmesh(itert)+x-xmesh(iterx), @(x) df(x)*tmesh(itert)+1 , rand, 10^(-4), 500 );
%         Newton.Imax(itert,iterx)  = length(x0);
%         Newton.Error(itert,iterx) = ex0(end);
%        u(itert,iterx)            = f(x0(end));
%         % Matlab fzero
%         x0                        = fzero(@(x) f(x)*tmesh(itert)+x-xmesh(iterx),rand);
%         Fzero.Error(itert,iterx)  = f(x0)*tmesh(itert)+x0-xmesh(iterx);
%         u(itert,iterx)            = f(x0(end));
%     end
% end
% clearvars x0 ex0

%% solve pde        
for itert = 2 : LenT
    for iterx = 2 : N
        maxu            = max(abs(u(itert-1,:)));
        Fplus           = 1/2*(1/2*(u(itert-1,iterx))^2+1/2*(u(itert-1,iterx+1))^2)+maxu/2*(u(itert-1,iterx)-u(itert-1,iterx+1));
        Fminus          = 1/2*(1/2*(u(itert-1,iterx))^2+1/2*(u(itert-1,iterx-1))^2)+maxu/2*(u(itert-1,iterx-1)-u(itert-1,iterx));
        diff            = ax(iterx)*(u(itert-1,iterx+1)-2*u(itert-1,iterx)+u(itert-1,iterx-1))/deltax^2;
        u(itert,iterx)  = u(itert-1,iterx) + deltat*(diff - (Fplus-Fminus)/(deltax));
    end
end
%u(:,3:5:end)= NaN;

clearvars maxu Fplus Fminus diff

% plot solution
figure;
tclick = 10;
plot(xmesh,u(1,:));hold on;title('Solution')
%saveas(gcf,'u.png')
for iter=1 : floor(LenT/tclick)
   plot(xmesh,u(iter*tclick,:)); hold on;
end
hold off
%saveas(gcf,'fig1.png')


% sample data
% DownFx         = 10;
% Dataxindex     = 1:DownFx:length(xmesh);
% Dataxmesh      = xmesh(Dataxindex);
% DataLenx       = length(Dataxmesh);
% Datadeltax     = Dataxmesh(2)-Dataxmesh(1);

i = 10;
gap = 10;
Dataxindex = [1];
j = 2;
while (i<length(xmesh))
    r = round(i + gap*rand(1));
    if r ~= Dataxindex(j-1)
        Dataxindex = [Dataxindex r];
        i = i+gap;
        j=j+1;
    end
end

Dataxmesh      = xmesh(Dataxindex);
DataLenx       = length(Dataxmesh);
Datadeltax     = Dataxmesh(2)-Dataxmesh(1);
        

DownFt         = 4;
Datatindex     = 1:DownFt:round(LenT);
Datatmesh      = tmesh(Datatindex);
DataT          = Datatmesh(end);
DataLent       = length(Datatindex);
Datadeltat     = Datatmesh(2)-Datatmesh(1);

Datau          = u(Datatindex,Dataxindex);

%% add gaussian noise
% sigma              = (noise/100)*norm(Datau,'fro')/sqrt(DataLent*DataLenx);
% Datau_noisy        = Datau + normrnd(0,sigma,DataLent,DataLenx);
% Datau_noisy(:,1)   = Bleft(0);
% Datau_noisy(:,end) = Bright(0);
% 
% 
% % plot noisy solution
% if noise>0
% figure;
% tclick = 4;
% plot(Dataxmesh,Datau_noisy(1,:)); hold on;title('Noisy solution')
% for iter=1 : floor(DataLent/tclick)
%    plot(Dataxmesh,Datau_noisy(iter*tclick,:)); hold on;
% end
% hold off
% saveas(gcf,'fig2_noisy.png')
% 
% end
Datau_noisy = Datau;
%% evaluate basis
LenB           = 1;
Phi            = ones(DataLenx,LenB);

% % 3-point stencil
% xrange = 2:DataLenx-1;
% dux    = ComputeDerivatives3P(Datau_noisy',Datadeltax);
% duxx   = ComputeDerivatives3P(dux,Datadeltax);

% %% 5-point stencil
% xrange = 3:DataLenx-2;
% dux    = ComputeDerivatives5P(Datau_noisy',Datadeltax);
% duxx   = ComputeDerivatives5P(dux,Datadeltax);
%  );
%% 5-point ENO
xrange = 5:DataLenx-4;
dux    = ComputeDerivativesEno5P(Datau',Datadeltax,Dataxmesh);
duxx   = ComputeDerivativesEno5P(dux,Datadeltax,Dataxmesh);




%% Form the linear system
LenD              = 10; % 10 dictionary atoms
uvec              = (Datau_noisy(2:end,xrange))';
uvec              = uvec(:);
dux               = dux(xrange,2:end);
dux               = dux(:);
duxx              = duxx(xrange,2:end);
duxx              = duxx(:);
A                 = [ones(length(dux),1) uvec  uvec.^2 ...
                      dux  dux.^2   uvec.*dux   duxx   duxx.^2     uvec.*duxx   dux.*duxx];  
dut               = ((Datau_noisy(2:end,xrange)-Datau_noisy(1:end-1,xrange))/Datadeltat)';
dut               = dut(:);
clear uvec dux duxx



% check the ground truth
alpha_true             = zeros(10,1);
supp                   = [6 7];
alpha_true(supp)      = [-1 couxx]';

% ResAbs = abs(dut-A*alpha_true);
% ResRe  = ResAbs./abs(dut);
% ResAbs = vec2mat(ResAbs,length(xrange));
% ResRe  = vec2mat(ResRe,length(xrange));
% 
% 
% figure;
% 
% %plot(1:length(du),abs(du-A*alpha_true),'red'); title('absolute error Ax-b')
% imagesc(ResAbs); colorbar;
% xlabel(['absolute residual of the linear system = ' num2str(norm(dut-A*alpha_true))])
% saveas(gcf,'fig3.png')

% % % % figure;
% % % % %plot(1:length(du),abs(du-A*alpha_true)./abs(du),'red'); title('relative error Ax-b')
% % % % imagesc(ResRe); colorbar;
% % % % xlabel(['relative residual of the linear system = ' num2str(norm(dut-A*alpha_true)/norm(du))])


% Acolnorm = sqrt(sum((abs(A)).^2,1));
% Anormed  = A*diag(1./Acolnorm);
% 
% 
% % look at coherence patter
% CohMatrix = abs(Anormed'*Anormed);
% figure
% imagesc(CohMatrix);colorbar
% title('Coherence pattern')
% saveas(gcf,'fig4_coherence.png')


b                   = dut;
A0                  = A;
Acolmax             = max(abs(A),[],1);
Amaxnorm            = A*diag(1./Acolmax);
alpha_truemaxnorm   = alpha_true.*Acolmax';

if IfNormCol == 0       % column not normalized
    A  = A0;
    alpha_trueshow = alpha_true;
elseif IfNormCol == 1     % normorlized by each column max
    A  = Amaxnorm;
    alpha_trueshow = alpha_truemaxnorm;
elseif IfNormCol == 2   % normlized by column norms
    A               = Anormed;
    alpha_trueshow  = alpha_true.*Acolnorm';
end

[RA , CA]           = size(A);

%% noise to signal coefficient with the min coefficient
conoise            = norm(b-A0*alpha_true)/min(abs(nonzeros(alpha_trueshow)));
fprintf('noise = %6.4f p ls-noise =%6.0f \n',noise,conoise)

% ADMM for group lasso 1/2||Ax-b||^2+lambda||x||_Group1
opts     = struct('itermax',200,'abstol',1e-4,'reltol',1e-2,'QUIET',1);
rho      = 1;

[xk , history] = ADMMGroupLasso(A , b , [], LenB , lambda, rho, opts);



[obj_true , objloss_true , objl1_true]  = objective(A, b, lambda, LenB , alpha_true);
fprintf('obj_true = %10.2f objloss_true = %10.2f objreg_true = %10.2f \n', obj_true,objloss_true,objl1_true);
[obj_rec , objloss_rec , objl1_rec]  = objective(A, b, lambda, LenB , xk);
fprintf('obj_rec  = %10.2f objloss_rec  = %10.2f objreg_rec  = %10.2f \n', obj_rec,objloss_rec,objl1_rec);



figure
subplot(1,2,1)
stem(1:CA,alpha_trueshow,'ro');hold on;
stem(1:CA,xk,'b*')
title('exact and recovered coefficients by ADMM')
subplot(1,2,2)
stem(1:CA,alpha_trueshow-xk,'ro');hold on;
title('difference')
saveas(gcf,'fig5.png')


% % % %%%%%%%%% norm x back
if IfNormCol == 1
    xk_normback = xk./Acolmax';
elseif IfNormCol == 2
    xk_normback = xk./Acolnorm';
end
figure
stem(1:CA,alpha_true,'ro');hold on;
stem(1:CA,xk_normback,'b*')
title('exact and recovered coefficients by ADMM')
legend('exact','recovered')
saveas(gcf,'fig6.png')

% subplot(1,2,2)
% stem(1:CA,alpha_true-xk_normback,'ro');hold on;
% title('difference')



