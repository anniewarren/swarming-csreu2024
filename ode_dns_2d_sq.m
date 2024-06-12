clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup (size of domain is 2pi x 2pi)
rng(1491283) % fix seed of rng to reproduce results

n = 4; % n = sqrt(N), where N is the total # of particles
N = n^2;
x = linspace(0, 2*pi, n+1); x = x(1:end-1)'; % erase last term

% X == X(j,k,1), Y == X(j,k,2), a.k.a all 2D positions
X = repmat(x,1,n);
Y = X';

% column by column into 1 vec for each -> p1 = (X(1), Y(1)) in 2D plane
% both length N = n^2
X = reshape(X,[],1);
Y = reshape(Y,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perturbations, parameters
dx = (X(2) - X(1))/n; % a step size in x, = (2*pi/n)/n = 2*pi/N
dy = (Y(n+1) - Y(1))/n; % a step size in y, should be 2*pi/N

epsilon = 0.3;  % width of potential

X = X + epsilon*cos(X);  % IC
X = X + 1e-2/N * randn(N,1); % adding some noise
%X(1:n) = mod(X(1:n) + .2, 2*pi);  % shift 1st row slightly

Y = Y + epsilon*cos(Y);    % IC
Y = Y + 1e-2/N * randn(N,1); % adding some noise

% true kappa = 1/pi^2?
%kappa = 0.0448;    
kappa = 0.053;   % 0.07 square vacuum?, n=8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting 

h=figure(1); % current figure

% animation speed
dt=.05;  %.1

j=0;
tplot=5;    % 50
nplot=round(tplot/dt); % only plot every nplot-th frame
%nplot=1;

% total time
t=0;
tend=4000; %30000  

while t<tend
    j=j+1; t=t+dt;

    F0 = V1_0(X, Y, N, epsilon, kappa);  % size Nx2
    X1 = X + dt*F0(:,1);    
    Y1 = Y + dt*F0(:,2);    
    F1 = V1_0(X1, Y1, N, epsilon, kappa); 
    VF = norm(F1,'inf')/dt; % current motion

    X = X + dt/2*(F0(:,1) + F1(:,1));
    Y = Y + dt/2*(F0(:,2) + F1(:,2));

    X = mod(X/(2*pi),1)*2*pi;
    Y = mod(Y/(2*pi),1)*2*pi;

    if mod(j,nplot)==0  
            set(0, 'currentfigure', h) % update the figure
            plot(X,Y,'.', MarkerSize=12)
            hold off
            
            title(['time t=' num2str(t) ', N=' num2str(N) ', VF=' num2str(VF) ', kappa=' num2str(kappa)])
            axis([ 0 2*pi 0 2*pi])
            xlabel('x')
            ylabel('y')
            drawnow
    end
    
    % stop animation once system settles
    if VF < 8e-4
        tend = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance matrices
% returns NxN matrices
function [Dx, Dy] = Dist(X,Y,N) 
    
    % x component distances mod 2pi
    Dx = mod(ones(N,1)*X' - X*ones(1,N), 2*pi);
    Dx = Dx - 2*pi*(Dx > pi);
    
    % y component distances mod 2pi
    Dy = mod(ones(N,1)*Y' - Y*ones(1,N), 2*pi);
    Dy = Dy - 2*pi*(Dy > pi);
    
    % euclidean distance
    %D = sqrt(Dx.^2 + Dy.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potential function
function G = V1_0(X, Y, N, epsilon, kappa)
    [Dx, Dy] = Dist(X,Y,N);
    G = zeros(N,2); % grad V at each point to be returned

    % Besselk function with 8 translates (modified bessel of 2nd kind)
    BK = (-1/(epsilon^3 * 2 * pi)) *...
        (...
              besselk(1, sqrt(Dx.^2 + Dy.^2)/epsilon)./sqrt(Dx.^2 + Dy.^2)...
            + besselk(1, sqrt((Dx+2*pi).^2 + Dy.^2)/epsilon)./sqrt((Dx+2*pi).^2 + Dy.^2)...
            + besselk(1, sqrt((Dx-2*pi).^2 + Dy.^2)/epsilon)./sqrt((Dx-2*pi).^2 + Dy.^2)...
            + besselk(1, sqrt(Dx.^2 + (Dy+2*pi).^2)/epsilon)./sqrt(Dx.^2 + (Dy+2*pi).^2)...
            + besselk(1, sqrt(Dx.^2 + (Dy-2*pi).^2)/epsilon)./sqrt(Dx.^2 + (Dy-2*pi).^2)...
            + besselk(1, sqrt((Dx+2*pi).^2 + (Dy+2*pi).^2)/epsilon)./sqrt((Dx+2*pi).^2 + (Dy+2*pi).^2)...
            + besselk(1, sqrt((Dx-2*pi).^2 + (Dy+2*pi).^2)/epsilon)./sqrt((Dx-2*pi).^2 + (Dy+2*pi).^2)...
            + besselk(1, sqrt((Dx-2*pi).^2 + (Dy-2*pi).^2)/epsilon)./sqrt((Dx-2*pi).^2 + (Dy-2*pi).^2)...
            + besselk(1, sqrt((Dx+2*pi).^2 + (Dy-2*pi).^2)/epsilon)./sqrt((Dx+2*pi).^2 + (Dy-2*pi).^2)...
        );
    
    % x components
    T = BK .* Dx + kappa * sin(Dx).*cos(Dy); % matrix
    T(~isfinite(T)) = 0; % fix dividing by 0 on diagonal of D, removes any NaN's
    G(:,1) = sum(T,2);  % column vector, k'th component is sum of k'th row of T

    % y components
    T = BK .* Dy + kappa * cos(Dx).*sin(Dy);
    T(~isfinite(T)) = 0;
    G(:,2) = sum(T,2);

    % average force acting on particle
    G = G * 1/N;
end