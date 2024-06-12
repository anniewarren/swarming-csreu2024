% direct simulation based on fft for swarming with repulsion (Dirac-like) and attraction (cos like)


close all;clear all;
%%%%%%%%%% numerical parameters %%%%%%%%%%%%%%%%%%%
dt = .003; % time stepping
n=256; % number of FOurier modes, best 2^N

eps=0.01;%artificial viscosity; small eps requires small dt
eps1=0.1; % smoothing of dirac conv kernel



%%%%%%%%%%%% system parameters %%%%%%%%%%%%%%%%%%
L = pi; % domain size
%xgrid = linspace(-L, L, n)';
xgrid = linspace(-L,L,n+1)';xgrid=xgrid(1:end-1);
mu=0.4;% bifurcation parameter, instability roughly at mu=1/pi~0.318 


%%%%%%%%%%%%%%%%%% convolution and derivative vectors%%%%%%%%%
k = [0:n/2 -(n/2-1):-1]';
k2 = k.^2; % not used here
% define various convolution kernels
cosconv=zeros(n,1);
cosconv(2)=1/2;
cosconv(end)=1/2; % this is the cosine kernel
delconv=ones(n,1);
delconv=1./(2*pi*(1+eps1.^2*k.^2)); %this is a smoothed version of the Dirac delta -- more modes if eps1 is small!
kon=1i*k.*(delconv -mu*cosconv); % putting the kernels together
visc = [1+dt*eps*k2]; % artificial viscosity

%%%%%%%%%%% %  initial shape %%%%%%%%%%%%%%%%%%%%%%%%
u0=1; %initial concentration
u = u0*ones(n,1)+0.1*(randn(n,1));u=u-(sum(u)/n-1);%random perturbation opreserving average
u = u0*ones(n,1)+0.7*sin(xgrid);%sine perturbation
uf = fft(u);


%%%%%%%%%  initialize time stepping %%%%%%%%%%%%%%%%%%%%%%%%%%
t_max=50000;
t_plot=1;
t=0;
h=figure(17);

while t < t_max
    t = t+dt;
    u = max(ifft(uf,'symmetric'),0); % Fourier might not preserve positivity, so hardcode this in here
    if mod(t,t_plot)<dt
        set(0,'CurrentFigure',h)
        plot(xgrid,u);
        axis([xgrid(1) xgrid(floor(end))  -0.1 3]);
        title(['Swarming, t = ' num2str(t) ', max(u)=' num2str(max(u))]);
        drawnow;
    end
    % now the actual time stepping, all the work done here
    uf =  (uf + dt*1i*k.*fft(u.*ifft(kon.*uf,'symmetric')))./visc;
end


 
