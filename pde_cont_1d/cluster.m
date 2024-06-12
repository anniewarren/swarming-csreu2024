clear all;
gpuon=0;
verbose=1;

%% problem parameters 
n=64*1024;% discretization size
M =1; % size of domain in units of 2pi
L = M*pi; 
x = linspace(-L, L, n+1)'; x=x(1:end-1); %n points total, end is periodic point
dx=2*L/n; 
%derivative vectors (in Fourier space)
k = ([[0:n/2] [-n/2+1: -1]]./M)';
k2=k.^2;
%  k4 = k2.^2;
Dlin = -k2 - 1;% linearization of potential differential operator (-1 prevents div by 0) 
pc= @(du) [ifft(fft(du(1:end-3))./Dlin,'symmetric');du(end-2:end)]; % preconditioner
D2=-k2;
D1=1i*k;

%initial shape
u0 = 1+7e-1*cos(x);
s=0; % s for drift (parameter compensating for extra phase condition)
alpha=0; % alpha for mass loss (mass constraint)
mu=-2.1; % mu bifurcation parameter
u=[u0;s;alpha;mu]; 

tang=[cos(x);0;0;0.01];tang=tang/norm(tang);

uinit=u;
% append parameters, ...
% s for drift, 
% alpha for mass loss, 
% mu for continuation parameter in secant

% initial regularizing viscosity
eps=0.001; % smaller epsilon requires more fourier modes
% potential
V1=-sin(x); % first derivative of only cos part of potential
V1h=2*pi*fft(V1)./n; % averaged Fourier Transform, so n is not needed in multiplication

% define objective function without parameters
spike= @(u) spikea(u,uinit,tang,D1,D2,eps,V1h);
% derivative will be called as 
% dspikea= @(du,u,uinit,tang,D1,D2,eps,V1h)

% Newton parameters
ntol=1e-4;
nmaxiter = 50;
nniter = 1;
nnitermax= 35;
nnincr=ones(n+1,1);
nminstep=-1e-8;
nrhs=0;
% initial res
nresidual=norm(spike(u));
ngmrestol=1e-6;


% arclength continuation initial stepsize
ds=1e-1;
contsteps=1000;
% compute residual
nrhs=spike(u);                     % compute residual
nresidual=norm(nrhs);           % estimate residual
du=u;

tic
% 
% if gpuon 
%     u=gpuArray(u);
% %      du=gpuArray(du);
%     x=gpuArray(x);
%     Dlin=gpuArray(Dlin);
%     nrhs=gpuArray(nrhs);
%     ngmrestol=gpuArray(ngmrestol);
% %      pc=gpuArray(pc);
%     nnincr=gpuArray(nnincr);
%     nniter=gpuArray(nniter);
%     nnitermax=gpuArray(nnitermax);
%     
%     sprintf('calculating in double precision on GPU')
% end
uinit=u;
% Newton loop
[u,newton_flag]=newtonloop(nresidual,ntol,du,u,uinit,tang,D1,D2,eps,V1h,nrhs,ngmrestol,pc,nnitermax,nminstep);

display(['first solution computed, mu=' num2str(u(end)) ', norm=' num2str(norm(u(1:end-3),'inf'))])
uold=u;
u=u-ds*tang;
% compute residual
nrhs=spike(u);                     % compute residual
nresidual=norm(nrhs);           % estimate residual

uinit=u;
[u,newton_flag]=newtonloop(nresidual,ntol,du,u,uinit,tang,D1,D2,eps,V1h,nrhs,ngmrestol,pc,nnitermax,nminstep);


display(['second solution computed, mu=' num2str(u(end)) ', norm=' num2str(norm(u(1:end-3),'inf'))])
tang=u-uold;tang=tang/norm(tang);

mu_list=[];
inf_list=[];
h=figure(1);
for contcount=1:contsteps
    uold=u;
    uinit=uold+ds*tang;
    u=uinit;
    spike= @(u) spikea(u,uinit,tang,D1,D2,eps,V1h);
    nrhs=spike(u);                     % compute residual
    nresidual=norm(nrhs);           % estimate residual
    [u,newton_flag]=newtonloop(nresidual,ntol,du,u,uinit,tang,D1,D2,eps,V1h,nrhs,ngmrestol,pc,nnitermax,nminstep);
    if newton_flag.error>0
        u=uold;
        ds=ds/2;
        display(['newton did not converge, reducing stepsize to ' num2str(ds)])
    else
        if verbose
            display([num2str(newton_flag.nniter) ' newton iterations'])
            display([num2str(newton_flag.nresidual) ' residual after newton'])
        end
        tang=u-uold;tang=tang/norm(tang);
        set(0,'CurrentFigure',h)
            plot(x,u(1:end-3));
            title(['$\kappa=$' num2str(u(end)) ', $|u|_\infty$=' num2str(norm(u(1:end-3),'inf')) ' $ds=$' num2str(ds)],'interpreter','latex')
            drawnow
        mu_list=[mu_list;u(end)];
        inf_list=[inf_list;norm(u(1:end-3),'inf')];
        if newton_flag.nniter >5
            ds=ds/2;
        end
        if newton_flag.nniter <3
            ds=min(1,ds*1.2);
        end
    end
    if max(u(1:end-3))>2.3
        break
    end
end
toc
figure(33)
plot(mu_list,inf_list,'.-')
xlabel(['$\mu$'],'Interpreter','latex')



% %  check Jacobian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  duh=rand(n,1)*1e-5;
%  res=RSH(uh+duh)-RSH(uh)-A(duh,uh,par);
%  plot(abs(res));
% %  gmres test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  sol=gmres(dspike,rand(n,1),10,1e-6,15,pc);%gmres(A,b,restart,tol,maxit,M) 
 

 
