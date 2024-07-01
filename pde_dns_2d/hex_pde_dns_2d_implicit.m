% direct simulation based on fft for 2D swarming with repulsion (Dirac delta) and attraction (cosine)

close all;
%%%%%%%%%% numerical parameters %%%%%%%%%%%%%%%%%%%
dt = 0.01; % time step
n = 128; % number of Fourier modes, best 2^N
iterations = 5000000;
plotintervals = 5*dt;
threshold = 0.1;
mu = 1/(sqrt(3)*pi*pi)-0.0075; % bifurcation parameter, instability roughly at 


%%%%%%%%%%%% system parameters %%%%%%%%%%%%%%%%%%
L = 2*pi; % domain size
y1grid = linspace(-L,L,n+1)';y1grid=y1grid(1:end-1); %column vector in hex
y2grid = linspace(-L,L,n+1);y2grid=y2grid(1:end-1); %row vector in hex
[Y1,Y2] = meshgrid(y1grid,y2grid); % mesh in hex


%%%%%%%%%%%%% change of basis components %%%%%%%%%%%%%
% hex basis vectors
e1=[1;0];
e2=[-1/2;sqrt(3)/2];
T = [e1,e2]; % square to hex change of basis matrix

% dual basis vectors
e1s = [1;1/sqrt(3)]; 
e2s = [0;2/sqrt(3)];
S = [e1s'; e2s']; % hex to square change of basis matrix

% calculating boundaries for plotting
lowerLeftbound = T*[-L;-L];
lowerRightbound = T*[L;-L];
upperLeftbound = T*[-L;L];
upperRightbound = T*[L;L];
ylowlim = lowerRightbound(2);
yuplim = upperRightbound(2);

% axes for plotting hex
X1 = T(1,1)*Y1 + T(1,2)*Y2;
X2 = T(2,1)*Y1 + T(2,2)*Y2; % here T(2,1) = 0

%%%%%%%%%%%%%%%%%% convolution and derivative vectors%%%%%%%%%
k = repmat([0:n/2 -n/2+1:-1]',1,n); % x fourier indices as a matrix
l = repmat([0:n/2 -n/2+1:-1],n,1); % y fourier indices as a matrix
% neglaplace = k.^2 + l.^2;
neglaplace = (4/3)*(k.^2 + k.*l + l.^2);
dy1 = 1i*k;
dy2 = 1i*l;
% define fourier coefficients of potential in square
vphys = cos(Y1)+cos(Y2)+cos(Y1-Y2);
v = fft2(vphys);

% v=zeros(n,n); 
% v(2,2)=(n*n)/4;
% v(2,n)=(n*n)/4;
% v(n,2)=(n*n)/4;
% v(n,n)=(n*n)/4; % cos(x)cos(y) fourier coefficients
    
%%%%%%%%%%% %  initial shape %%%%%%%%%%%%%%%%%%%%%%%%
u0=1; %initial concentration
%u = u0*ones(n,n)+0.1*(randn(n,n));
%u=u.*((n*n)/sum(sum(u))); %random perturbation preserving average
%u = u0*ones(n,n)+0.7*sin(xgrid)*sin(ygrid); %sine perturbation

u = u0*ones(n,n)+0.1*(cos(Y1)+cos(Y2)+cos(Y1-Y2));
u=u.*((n*n)/sum(sum(u))); %preserving average
uf = fft2(u);

%%%%%%%%%  initialize time stepping %%%%%%%%%%%%%%%%%%%%%%%%%%
t=0;
t_max = plotintervals*iterations;
h=figure(Name='Simulation');
iter = 1;

while t < t_max
    maxiter= 50;
    t = t+dt;
    u = max(ifft2(uf,'symmetric'),0); % force positivity of inv Fourier transform
    
    if mod(t,plotintervals)<dt % only plot at intervals of t_plot
        set(0,'CurrentFigure',h)
        parulaedited = [1 1 1;parula];
        

        subplot(2,1,1)
        surf(X1,X2,u,EdgeColor="none"); %3D surface plot
        view(2);
        title(['Swarming, t = ' num2str(t) ', max(u)=' num2str(max(max(u))), ', mass=' num2str(sum(sum(u)))]);
        axis([upperLeftbound(1) lowerRightbound(1) ylowlim yuplim 0.1 6])
        daspect([1 1 1])
        clim([0.1,6]);
        colorbar;

        subplot(2,1,2)
        surf(X1,X2,u,EdgeColor="none"); %3D surface plot
        axis([upperLeftbound(1) lowerRightbound(1) ylowlim yuplim 0.1 6])
        daspect([1 1 1])
        clim([0.1,6]);
        colorbar;
        drawnow;
        %saveas(gcf,['swarmpics/n' num2str(n) '_mu' sprintf('%.4f', mu) '_iter' sprintf('%04d',iter) '.jpg']);
        iter = iter + 1;
    end
    
    % now the actual time stepping, all the work done here
    xcom = fft2(u.*ifft2(dy1.*(v.*uf.*(4*pi*pi/(n*n))),'symmetric'));
    ycom = fft2(u.*ifft2(((1/sqrt(3))*dy1 + (2/sqrt(3))*dy2).*(v.*uf.*(4*pi*pi/(n*n))),'symmetric'));
    F = fft2(-1*mu*ifft2(dy1.*xcom + ((1/sqrt(3))*dy1 + (2/sqrt(3))*dy2).*ycom,'symmetric'));
    Fvec = reshape(F,n*n,1);

    D0 =@(w) w - dt.*D(w,uf,n,dy1,dy2); %takes in a vector w, returns a vector
    
    M =@(v) spdiags(reshape((ones(n,n) + dt.*neglaplace).^(-1),n*n,1),0,n*n,n*n)*v; % precondition matrix
    
    ufvec = reshape(uf,n*n,1);
    disp(['time ' num2str(t)])
    ufvecnew = cgs(D0,ufvec+dt.*Fvec,[],maxiter,M);
    % while convergence == 1 %check that cgs converged, if not try increasing iterations
    %     if maxiter < 100
    %         maxiter = 100;
    %         [ufvecnew,convergence] = cgs(D0,ufvec+dt.*Fvec,[],maxiter,M);
    %     else
    %         return %if increasing iterations doesn't help, break
    %     end
    % end
    %norm(ufvecnew-ufvec);
    uf = reshape(ufvecnew,n,n); % finally convert back to a matrix
end

% takes in un as an nxn matrix, un1 an n^2 column vector both in fourier space 
% outputs an n^2 column vector in fourier space
function f = D(un1vec,un,n,dy1,dy2) 
    epsilon=1e-6; % artificial viscosity; small eps requires small dt
    % dy1 = dy1;
    % dy2 = ((1/sqrt(3))*dy1 + (2/sqrt(3))*dy2);
    un1 = reshape(un1vec,n,n);
    y1com = fft2((ifft2(un,'symmetric')+epsilon.*ones(n,n)).*ifft2(dy1.*un1, 'symmetric'));
    y2com = fft2((ifft2(un,'symmetric')+epsilon.*ones(n,n)).*ifft2(((1/sqrt(3))*dy1 + (2/sqrt(3))*dy2).*un1, 'symmetric'));
    fmatrix = dy1.*y1com + ((1/sqrt(3))*dy1 + (2/sqrt(3))*dy2).*y2com;
    f = reshape(fmatrix,n*n,1);
end