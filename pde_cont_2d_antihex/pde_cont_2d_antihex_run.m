function [mu_list_out,vac_list_out] = pde_cont_2d_antihex_run(nFouriermodes,epsilon,verbose,plotting,gpuON)
    %verbose = true;
    %plotting = true;
    %gpuON = false;
    
    % set up parameters
    n=nFouriermodes;% discretization size
    eps=epsilon; % regularizing viscosity- smaller eps requires larger n ~ 1/sqrt(eps)
    threshold = eps;
    ds=1e-3; % secant step size
    M =2; % size of domain in units of 2pi
    L = M*pi; 
    y1 = linspace(-L, L, n+1)'; y1=y1(1:end-1); %n points total, end is periodic point
    y2 = linspace(-L,L,n+1);y2=y2(1:end-1);
    %dx=2*L/n; 
    %dy=2*L/n;
    
    [Y1,Y2]= meshgrid(y1,y2);
    
    %%% hexagonal lattice info
    e1=[1;0];
    e2=[-1/2;sqrt(3)/2];
    % crazy train out of crazy town %
    T = [e1,e2];
    
    e1s = [1;1/sqrt(3)];
    e2s = [0;2/sqrt(3)];
    
    % calculating boundaries
    lowerLeftbound = T*[-L;-L];
    lowerRightbound = T*[L;-L];
    upperLeftbound = T*[-L;L];
    upperRightbound = T*[L;L];
    ylowlim = lowerRightbound(2);
    yuplim = upperRightbound(2); 
    
    % project onto square axes for plotting
    X1 = T(1,1)*Y1 + T(1,2)*Y2;
    X2 = T(2,1)*Y1 + T(2,2)*Y2; % here T(2,1) = 0
    
    %derivative vectors (in Fourier space)
    k = repmat([0:n/2 -n/2+1:-1]',1,n)./M; % x fourier indices as a matrix
    l = repmat([0:n/2 -n/2+1:-1],n,1)./M; % y fourier indices as a matrix
    D1x1=1i*k;
    D1x2=1i*l;
    D2x1 = -k.^2;
    D2x2 = -l.^2;
    laplace = (4/3)*(-k.^2 - k.*l - l.^2);
    
    % potential
    V = cos(Y1)+cos(Y2)+cos(Y1-Y2);
    Vhat = fft2(V).*((2*sqrt(3)*pi*pi)/(n*n));
    
    %initial conditions
    u0=1; %initial concentration
    uinit = u0*ones(n,n)-0.1*(cos(Y1)+cos(Y2)+cos(Y1-Y2));
    uinit=uinit.*((n*n)/sum(sum(uinit))); %preserving average
    
    %%% GPU data transfers %%%
    if gpuON
        n = gpuArray(n);
        eps = gpuArray(eps);
        threshold = gpuArray(threshold);
        D1x1 = gpuArray(D1x1);
        D1x2 = gpuArray(D1x2);
        laplace = gpuArray(laplace);
        Vhat = gpuArray(Vhat);
        uvec = gpuArray(uvec);
        tangentvec = gpuArray(tangentvec);
        ntol = gpuArray(ntol);
        ngmrestol = gpuArray(ngmrestol);
        nnitermax = gpuArray(nnitermax);
        nminstep = gpuArray(nminstep);
        cutoff = gpuArray(cutoff);
        contsteps = gpuArray(contsteps);
        errorsinarow = zeros(1,1,"gpuArray");
        mu_list = zeros(1,contsteps,"gpuArray");
        vac_list = zeros(1,contsteps,"gpuArray");
    end
    
    % preconditioner for gmres
    %Dlin = spdiags(reshape((laplace-1),n*n,1),0,n*n,n*n); % linearization of potential differential operator (-1 prevents div by 0) 
    pc= @(duvec) [reshape(ifft2(fft2(reshape(duvec(1:end-4),n,n))./(laplace-1),'symmetric'),n*n,1);duvec(end-3:end)]; % preconditioner
    
    
    if plotting
        figure(11)
        surf(X1,X2,uinit,EdgeColor="none")
            view(2)
            axis([upperLeftbound(1) lowerRightbound(1) ylowlim yuplim threshold 5])
            daspect([1 1 1])
            drawnow
    end
    
    mu=1/(sqrt(3)*pi*pi); % mu bifurcation parameter
    sx=0; % s for drift (compensating for extra phase condition), in x and y
    sy=0; 
    alpha=0; % alpha for mass loss 
    
    uvec = [reshape(uinit,n*n,1);sx;sy;alpha;mu];
    
    tangextraparams = [0 0 0 0.01];
    tangent=cos(Y1)+cos(Y2)+cos(Y1-Y2);
    tangentvec=[reshape(tangent,n*n,1);tangextraparams'];tangentvec=tangentvec/norm(tangentvec);
    
    % Newton parameters
    ntol=1e-4;
    ngmrestol=1e-6;
    nnitermax= 35;
    nminstep=-1e-8;
    
    %%%% GPU commands %%%%
    
    % do first newton
    duvec=uvec;
    u0vec=uvec;
    [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x1,D1x2,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
    disp(['first solution computed, mu=' num2str(uvec(end)) ', norm=' num2str(norm(uvec(1:end-4),'inf'))])
    uoldvec=uvec; % save first solution
    
    % do second newton
    uvec=uvec-ds*tangentvec; % take step
    u0vec=uvec;
    [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x1,D1x2,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
    disp(['second solution computed, mu=' num2str(uvec(end)) ', norm=' num2str(norm(uvec(1:end-4),'inf'))])
    tangentvec=uvec-uoldvec;tangentvec=tangentvec/norm(tangentvec);
    
    mu_list=[];
    %inf_list=[];
    vac_list=[];
    if plotting
        h=figure(1);
    end    
    %w=waitbar(0,'forming vacuum...');
    
    contsteps = 100000;
    vacuum_size = 0;
    
    %start secant continuation
    fprintf(['computing iteration... mu=' num2str(uvec(end),'%06.5f') ' umax=' num2str(norm(uvec(1:end-4),'inf'),'%06.5f') ' (vac: ' num2str(vacuum_size,'%06.3f') ')'])
    for contcount=1:contsteps
        uoldvec=uvec;
        u0vec=uoldvec+ds*tangentvec;
        uvec = u0vec;
        F= @(uvec) funk(uvec,u0vec,tangentvec,Vhat,D1x1,D1x2,laplace,eps,n);
        [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x1,D1x2,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
        
        if newtonflag.error>0
            uvec=uoldvec;
            ds=ds/2;
            display(['newton did not converge, reducing stepsize to ' num2str(ds)])
        else
            if verbose
                disp([num2str(newtonflag.iter) ' newton iterations'])
                disp([num2str(newtonflag.nresidual) ' residual after newton'])
            end
            tangentvec=uvec-uoldvec;tangentvec=tangentvec/norm(tangentvec);
    
            if plotting
                set(0,'CurrentFigure',h)
                u = reshape(uvec(1:end-4),n,n);
                subplot(2,1,1)
                surf(X1,X2,u,EdgeColor="none")
                view(2)
                axis([upperLeftbound(1) lowerRightbound(1) ylowlim yuplim threshold 5])
                daspect([1 1 1])
                title(['$\mu=$' num2str(uvec(end)) ', $|u|_\infty$=' num2str(norm(uvec(1:end-4),'inf')) ' $ds=$' num2str(ds)],'interpreter','latex')
                subplot(2,1,2)
                surf(X1,X2,u,EdgeColor="none")
                axis([upperLeftbound(1) lowerRightbound(1) ylowlim yuplim threshold 5])
                daspect([1 1 1])
                drawnow
            end    
    
            mu_list=[mu_list;uvec(end)];
            %inf_list=[inf_list;norm(uvec(1:end-4),'inf')];
    
            % calculate size of vacuum
            vacuum = uvec(1:end-4) < threshold;
            vacuum_size = sum(vacuum)*4*pi*2*sqrt(3)*pi/(n*n);
            if vacuum_size > 0
                %vac_pts = vac_pts+1;
                %waitbar(vac_pts/vac_pts_max,w,'taking vacuum data...')
            end
            %vac_radius = sqrt(vacuum_size/(2*pi));
            %disp(['vacuum size...' num2str(vacuum_size)])
            vac_list = [vac_list;vacuum_size];
            %l_list = [l_list;vac_radius];
            fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' num2str(uvec(end),'%06.5f') ' umax=' num2str(norm(uvec(1:end-4),'inf'),'%06.5f') ' (vac: ' num2str(vacuum_size,'%06.3f') ')'])
            if newtonflag.iter >5
                ds=ds/2; % decrease step size if newton not converging
            end
            if newtonflag.iter <3
                ds=min(1,ds*1.2); % increase stepsize if newton performing well
            end
        end
        if (vacuum_size/8)^2 > 0.8
            disp("reached breakpoint")
            break
        end
    end
    if gpuON
        mu_list_out = gather(mu_list);
        vac_list_out = gather(vac_list);
    else
        mu_list_out = mu_list;
        vac_list_out = vac_list;
    end
end