function [mu_list_out, vac_list_out, inf_list_out] = cpu_pde_cont_2d(nfouriermodes,epsilon)
    % set up parameters
    n=nfouriermodes;% discretization size
    eps=epsilon; % regularizing viscosity- smaller eps requires larger n ~ 1/sqrt(eps)
    verbose = false; 
    threshold = eps*2; %threshold for plotting and calculating vacuum size 
    ds=1e-1; % secant step size
    M =1; % size of domain in units of 2pi
    L = M*pi; 
    x = linspace(-L, L, n+1)'; x=x(1:end-1); %n points total, end is periodic point
    y = linspace(-L,L,n+1);y=y(1:end-1);
    %dx=2*L/n; 
    %dy=2*L/n;
    
    % arguments for displaying plots and saving them
    plotting = true;
    plotinterval = 20;
    gettingVideo = true;
    
    %derivative vectors (in Fourier space)
    k = repmat([0:n/2 -n/2+1:-1]',1,n)./M; % x fourier indices as a matrix
    l = repmat([0:n/2 -n/2+1:-1],n,1)./M; % y fourier indices as a matrix
    D1x=1i*k;
    D1y=1i*l;
    D2x = -k.^2;
    D2y = -l.^2;
    laplace = -k.^2 -l.^2;
    
    % potential
    V = cos(x)*cos(y);
    Vhat = fft2(V).*((4*pi*pi)/(n*n));
    
    % preconditioner for gmres
    pc= @(duvec) [reshape(ifft2(fft2(reshape(duvec(1:end-4),n,n))./(laplace-1),'symmetric'),n*n,1);duvec(end-3:end)]; % preconditioner
    
    %initial conditions
    %uinit = 1+7e-1*cos(x)*cos(y); % initial shape (past bifurcation)
    uinit = ones(n,n); % (before bifurcation)
    mu=0.1013-epsilon; % mu bifurcation parameter
    sx=0; % s for drift (compensating for extra phase condition), in x and y
    sy=0; 
    alpha=0; % alpha for mass loss 
    
    uvec = [reshape(uinit,n*n,1);sx;sy;alpha;mu];
    
    tangextraparams = [0 0 0 0.01];
    tangent=cos(x)*cos(y);
    tangentvec=[reshape(tangent,n*n,1);tangextraparams'];tangentvec=tangentvec/norm(tangentvec);
    
    % Newton parameters
    ntol=1e-4;
    ngmrestol=1e-6;
    nnitermax= 35;
    nminstep=-1e-8;
    
    
    
    % do first newton
    duvec=uvec;
    u0vec=uvec;
    [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x,D1y,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
    disp(['first solution computed, mu=' num2str(uvec(end)) ', norm=' num2str(norm(uvec(1:end-4),'inf'))])
    uoldvec=uvec; % save first solution
    
    % do second newton
    uvec=uvec+0.1*ds*tangentvec; % take step, 0.1 to keep step below bif to resolve branching off from crystal state
    u0vec=uvec;
    [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x,D1y,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
    disp(['second solution computed, mu=' num2str(uvec(end)) ', norm=' num2str(norm(uvec(1:end-4),'inf'))])
    tangentvec=uvec-uoldvec;tangentvec=tangentvec/norm(tangentvec);
    
    mu_list=[];
    inf_list=[];
    vac_list=[];
    l_list=[];
    
    if plotting 
        h=figure(1);
    end    

    contsteps = 1000000;
    vacuum_size=0;
   
    iter = 1;

    %start secant continuation
    vactrue = 'f';
    fprintf(['computing iteration... mu=' num2str(uvec(end),'%06.5f') ' umax=' num2str(norm(uvec(1:end-4),'inf'),'%06.5f') ' (vac: ' num2str(vacuum_size,'%06.3f') ')'])
    for contcount=1:contsteps
        uoldvec=uvec;
        u0vec=uoldvec+ds*tangentvec;
        uvec = u0vec;
        F= @(uvec) funk(uvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n);
        [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x,D1y,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
        
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
                if mod(contcount,plotinterval) == 1
                    set(0,'CurrentFigure',h)
                    u = reshape(uvec(1:end-4),n,n);
                    surf(x,y,u,LineStyle="none");
                    view(2)
                    title(['$\mu=$' num2str(uvec(end)) ', $|u|_\infty$=' num2str(norm(uvec(1:end-4),'inf')) ' $ds=$' num2str(ds)],'interpreter','latex')
                    axis([-pi pi -pi pi threshold 6])
                    xlabel('x')
                    ylabel('y')
                    zlabel('u')
                    drawnow
                    if gettingVideo
                        saveas(gcf,['bif_diag_inset/n' num2str(n) '_mu' sprintf('%.4f', uvec(end)) '_iter' num2str(iter,'%03u') '.fig']);
                    end
                    if vacuum_size > 0
                        plotinterval = 5; % get more frequent pics once things get interesting
                    end    
                    iter = iter+1;
                end
            end    
            mu_list=[mu_list;uvec(end)];
            inf_list=[inf_list;norm(uvec(1:end-4),'inf')];
    
            % calculate size of vacuum
            vacuum = uvec(1:end-4) < threshold;
            vacuum_size = sum(vacuum)*4*pi*pi/(n*n);
            if vacuum_size > 0
                vactrue = 't';
            end
            vac_radius = sqrt(vacuum_size/(2*pi));
            if verbose 
                disp(['vacuum size...' num2str(vacuum_size)])
            end    
            vac_list = [vac_list;vacuum_size];
            l_list = [l_list;vac_radius];
            
            fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' num2str(uvec(end),'%06.5f') ' umax=' num2str(norm(uvec(1:end-4),'inf'),'%06.5f') ' (vac: ' num2str(vacuum_size,'%06.3f') ')'])
            
            if newtonflag.iter >5
                ds=ds/2; % decrease step size if newton not converging
            end
            if newtonflag.iter <3
                ds=min(1,ds*1.2); % increase stepsize if newton performing well
            end
        end
        if vacuum_size > 5 % stopping point - can change to whatever is desired
            display("reached breakpoint")
            break
        end
    end
    inf_list_out = inf_list;
    l_list_out = l_list; % could pass this out if desired, but redundant given vac_list
    mu_list_out = mu_list;
    vac_list_out = vac_list;
end