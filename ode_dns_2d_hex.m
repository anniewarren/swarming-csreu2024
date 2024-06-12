clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
rng(81923756) % fix seed of rng to reproduce results

% hex basis vectors (e1, e2)
e1 = [1;0]; e2 = [-1/2; sqrt(3)/2];

% dual lattice vectors (e1*, e2*)
e1s = [1; 1/sqrt(3)]; e2s = [0; 2/sqrt(3)];

% n = sqrt(N), where N is the total # of particles
n = 5; 
N = n^2;

% X == W(j,k,1), Y == W(j,k,2), a.k.a all 2D positions
X = zeros(n); Y = zeros(n);
for i = 1:n
    for j = 1:n
        X(i,j) = (2*pi/n)*( i*e1(1) + j*e2(1) );
        Y(i,j) = (2*pi/n)*( i*e1(2) + j*e2(2) );
    end
end
% column by column into 1 vec for each -> p1 = (X(1), Y(1)) in 2D plane
% both length N = n^2
X = reshape(X,[],1);
Y = reshape(Y,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perturbations, parameters
dx = (X(2) - X(1))/n; % a step size in x, = 2*pi/N
dy = (Y(n+1) - Y(1))/n; % a step size in y, = sqrt(3)*pi/N

epsilon = 0.3;  % width of potential

%X = X - epsilon*cos(X);  % IC
X = X + 1e-2/N * randn(N,1); % adding some noise
%X(1:n) = mod(X(1:n) + .2, 2*pi);  % shift 1st row slightly

%Y = Y + epsilon*cos(Y);    % IC
Y = Y + 1e-2/N * randn(N,1); % adding some noise

% true kappa = 1/pi^2?
%kappa = 0.0448;    
kappa = 0.05; 

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

    % trapezoid rule
    F0 = V1_0(X, Y, N, epsilon, kappa, e1, e2, e1s, e2s);  % size Nx2
    X1 = X + dt*F0(:,1);    
    Y1 = Y + dt*F0(:,2);    
    F1 = V1_0(X1, Y1, N, epsilon, kappa, e1, e2, e1s, e2s); 
    VF = norm(F1,'inf')/dt; % current motion

    X = X + dt/2*(F0(:,1) + F1(:,1));
    Y = Y + dt/2*(F0(:,2) + F1(:,2));
    
    % project back into rhombus domain
    % [x y] = w*e1 + z*e2
    w = X + 1/sqrt(3) * Y;
    z = 2/sqrt(3) * Y;
    w = mod(w, 2*pi);
    z = mod(z, 2*pi);
    X = w - 1/2 * z;
    Y = sqrt(3)/2 * z;

    if mod(j,nplot)==0  
            set(0, 'currentfigure', h) % update the figure
            plot(X,Y,'.', MarkerSize=12)
            hold off
            
            title(['time t=' num2str(t) ', N=' num2str(N) ', VF=' num2str(VF) ', kappa=' num2str(kappa)])
            axis([-pi 2*pi 0 sqrt(3)*pi])
            xlabel('x')
            ylabel('y')
            drawnow
    end
    
    % stop animation once system settles
    if VF < 8e-4 && t > 1000
        tend = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance matrices
% returns NxN matrices
function [Dx, Dy] = Dist(X,Y,N) 
    
    % x,y component distances
    Dx = ones(N,1)*X' - X*ones(1,N);
    Dy = ones(N,1)*Y' - Y*ones(1,N);
    
    % re-center differences (bessel not really periodic)
    Dx = Dx + pi*(Dy>(sqrt(3)*pi/2)) - pi * (Dy < -sqrt(3)*pi/2);
    Dy = Dy - sqrt(3)*pi * (Dy > (sqrt(3)*pi/2)) + sqrt(3)*pi * (Dy < -sqrt(3)*pi/2);
    Dx = Dx - 2*pi*(Dx > pi) + 2*pi * (Dx < -pi);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potential function - different attractive part than square lattice
function G = V1_0(X, Y, N, epsilon, kappa, e1, e2, e1s, e2s)
    [Dx, Dy] = Dist(X,Y,N);
    G = zeros(N,2); % grad V at each point to be returned

    % Besselk function with 8 translates (modified bessel of 2nd kind)
    % (0,0), +-e1, +-e2, +-(e1+e2), +-(e1-e2)
    e1 = 2*pi * e1;
    e2 = 2*pi * e2;
    c1 = e1 + e2; 
    c2 = e1 - e2;

    BK = (-1/(epsilon^3 * 4 * pi^2)) *...
        (...
              besselk(1, sqrt(Dx.^2 + Dy.^2)/epsilon)./sqrt(Dx.^2 + Dy.^2)...
            + besselk(1, sqrt((Dx+e1(1)).^2 + (Dy+e1(2)).^2)/epsilon)./sqrt((Dx+e1(1)).^2 + (Dy+e1(2)).^2)...
            + besselk(1, sqrt((Dx-e1(1)).^2 + (Dy-e1(2)).^2)/epsilon)./sqrt((Dx-e1(1)).^2 + (Dy-e1(2)).^2)...
            + besselk(1, sqrt((Dx+e2(1)).^2 + (Dy+e2(2)).^2)/epsilon)./sqrt((Dx+e2(1)).^2 + (Dy+e2(2)).^2)...
            + besselk(1, sqrt((Dx-e2(1)).^2 + (Dy-e2(2)).^2)/epsilon)./sqrt((Dx-e2(1)).^2 + (Dy-e2(2)).^2)...
            + besselk(1, sqrt((Dx+c1(1)).^2 + (Dy+c1(2)).^2)/epsilon)./sqrt((Dx+c1(1)).^2 + (Dy+c1(2)).^2)...
            + besselk(1, sqrt((Dx-c1(1)).^2 + (Dy-c1(2)).^2)/epsilon)./sqrt((Dx-c1(1)).^2 + (Dy-c1(2)).^2)...
            + besselk(1, sqrt((Dx+c2(1)).^2 + (Dy+c2(2)).^2)/epsilon)./sqrt((Dx+c2(1)).^2 + (Dy+c2(2)).^2)...
            + besselk(1, sqrt((Dx-c2(1)).^2 + (Dy-c2(2)).^2)/epsilon)./sqrt((Dx-c2(1)).^2 + (Dy-c2(2)).^2)...
        );
    
    % x components
    T = BK .* Dx + kappa * (sin(Dx+e1s(2)*Dy) + sin(Dx + (e1s(2)-e2s(2))*Dy)); % matrix
    T(~isfinite(T)) = 0; % fix dividing by 0 on diagonal of D, removes any NaN's
    G(:,1) = sum(T,2);  % column vector, k'th component is sum of k'th row of T

    % y components
    T = BK .* Dy + kappa * (e1s(2)*sin(Dx+e1s(2)*Dy) + e2s(2)*sin(e2s(2)*Dy) + (e1s(2)-e2s(2))*sin(Dx + (e1s(2)-e2s(2))*Dy) );
    T(~isfinite(T)) = 0;
    G(:,2) = sum(T,2);

    % average force acting on particle
    G = G * 1/N;
end