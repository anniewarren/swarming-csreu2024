N = 131; % number of particles and 
%L=2*pi; % not in use but its size of domain
rng(1) % fix seed of rng to reproduce results
x = linspace(0, 2*pi, N+1); x = x(1:end-1)'; % erase last term
x0 = x';
dx = x0(2) - x0(1); % a step size in x. should be 2*pi/N
epsilon = 5*dx;
x = x + epsilon/10*cos(x); 
x = x + 1e-2/N*randn(N,1); % adding some noise
mu = 1.2112 * epsilon^2/pi;
% mu=.1;

V1_0 = @(y) -1/2*sign(y).*sinh((pi-abs(y))/epsilon)./sinh(pi/epsilon) + mu*sin(y); % potential function
F = @(x) sum(V1_0(Dist(x,N)), 2)/N; % average force acting on particle
% if A is a matrix, then sum(A,2) returns a column vector containing the sum of each row.

dt=.1;

t=0;
j=0;

tplot=50;
nplot=round(tplot/dt); % only plot every nplot-th frame
tend=30000;

% 
% nplot=1;

h=figure(1); % current figure

while t<tend
    j=j+1; t=t+dt;
    F0 = F(x);
    x1 = x + dt*F0;
    F1 = F(x1);
    x = x + dt/2*(F0+F1);
    x = mod(x/(2*pi),1)*2*pi;

    if mod(j,nplot)==0
            set(0, 'currentfigure', h) % update the figure
            axis([1 N 0 2*pi])
            plot(x,'.')
            hold off
            
            cg=sum(x)/N; % mean
            sd=sum((x-cg).^2)/N; % std
            title(['time t=' num2str(t) ', avgpos=' num2str(cg) ', sdeviation=' num2str(sd)])
            drawnow
    end
end

function F = Dist(x,N)
    % returns a matrix of size NxN whose entires a_ij are the distances
    % between partical i and partical j, mod 2*pi 
    % N must be the size of x
    F = mod(ones(N,1)*x' - x*ones(1,N), 2*pi);
    F = F - 2*pi*(F > pi);
end

