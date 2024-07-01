% direct simulation based on fft for 2D swarming with repulsion (Dirac delta) and attraction (cosine)
    
function pde_dns_implicit_amw(mu,nfouriermodes)
    close all;
    %%%%%%%%%% numerical parameters %%%%%%%%%%%%%%%%%%%
    dt = 0.001; % time step
    n = nfouriermodes; % number of Fourier modes, best 2^N
    iterations = 500;
    plotintervals = 0.01;
 
    
    %%%%%%%%%%%% system parameters %%%%%%%%%%%%%%%%%%
    L = pi; % domain size
    xgrid = linspace(-L, L, n+1)';xgrid=xgrid(1:end-1); %column vector
    ygrid = linspace(-L,L,n+1);ygrid=ygrid(1:end-1); %row vector
    mu; % bifurcation parameter, instability roughly at mu=1/pi^2 (0.1013)
    
    %%%%%%%%%%%%%%%%%% convolution and derivative vectors%%%%%%%%%
    k = repmat([0:n/2 -n/2+1:-1]',1,n); % x fourier indices as a matrix
    l = repmat([0:n/2 -n/2+1:-1],n,1); % y fourier indices as a matrix
    neglaplace = k.^2 + l.^2;
    dx = 1i*k;
    dy = 1i*l;
    % define fourier coefficients of potential
    vphys = repmat(cos(1.2*xgrid),1,n) + repmat(cos(1.2*ygrid),n,1);
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

    u = u0*ones(n,n)+0.9*repmat(sin(1.2*xgrid),1,n) + 0.5*repmat(sin(1.2*ygrid),n,1);
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
            subplot(2,1,1)
            surf(xgrid,ygrid,u,EdgeColor="none"); %3D surface plot
            subplot(2,1,2)
            surface(xgrid,ygrid,u,EdgeColor='none'); %2D surface plot
            axis([-pi pi -pi pi 0.001 2]);
            clim([0.001,2]);
            colorbar;
            title(['Swarming, t = ' num2str(t) ', max(u)=' num2str(max(max(u))), ', mass=' num2str(sum(sum(u)))]);
            drawnow;
            %saveas(gcf,['swarmpics/n' num2str(n) '_mu' sprintf('%.4f', mu) '_iter' sprintf('%04d',iter) '.jpg']);
            iter = iter + 1;
        end
        
        % now the actual time stepping, all the work done here
        xcom = fft2(u.*ifft2(dx.*(v.*uf.*(4*pi*pi/(n*n))),'symmetric'));
        ycom = fft2(u.*ifft2(dy.*(v.*uf.*(4*pi*pi/(n*n))),'symmetric'));
        F = fft2(-1*mu*ifft2(dx.*xcom + dy.*ycom,'symmetric')); %is this correct? changed mult by mu to phys space
        Fvec = reshape(F,n*n,1);
    
        D0 =@(w) w - dt.*D(w,uf,n,dx,dy); %takes in a vector w, returns a vector
        
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
end

% takes in un as an nxn matrix, un1 an n^2 column vector both in fourier space 
% outputs an n^2 column vector in fourier space
function f = D(un1vec,un,n,dx,dy) 
    epsilon=1e-3; % artificial viscosity; small eps requires small dt
    un1 = reshape(un1vec,n,n);
    xcom = fft2((ifft2(un,'symmetric')+epsilon.*ones(n,n)).*ifft2(dx.*un1, 'symmetric'));
    ycom = fft2((ifft2(un,'symmetric')+epsilon.*ones(n,n)).*ifft2(dy.*un1, 'symmetric'));
    fmatrix = dx.*xcom + dy.*ycom;
    f = reshape(fmatrix,n*n,1);
end