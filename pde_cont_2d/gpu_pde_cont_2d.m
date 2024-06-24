function [mu_list_out, vac_list_out] = gpu_pde_cont_2d(nfouriermodes,epsilon,verbose,contsteps)
    % set up parameters
    n=nfouriermodes;% discretization size
    eps=epsilon; % regularizing viscosity- smaller eps requires larger n ~ 1/sqrt(eps)
    ds=1e-1; % secant step size
    M =1; % size of domain in units of 2pi
    L = M*pi; 
    x = linspace(-L, L, n+1)'; x=x(1:end-1); %n points total, end is periodic point
    y = linspace(-L,L,n+1);y=y(1:end-1);
    %dx=2*L/n; 
    %dy=2*L/n;
    
    %derivative vectors (in Fourier space)
    k = repmat([0:n/2 -n/2+1:-1]',1,n)./M; % x fourier indices as a matrix
    l = repmat([0:n/2 -n/2+1:-1],n,1)./M; % y fourier indices as a matrix
    D1x=1i*k;
    D1y=1i*l;
    D2x = -k.^2;
    D2y = -l.^2;
    laplace = -1*(k.^2+l.^2);
    
    % potential
    V = cos(x)*cos(y);
    Vhat = fft2(V).*((4*pi*pi)/(n*n));
    
    %initial conditions
    uinit = 1+9e-1*cos(x)*cos(y); % initial shape
    mu=0.1013+eps; % mu bifurcation parameter
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
    
    % extra variables needed later
    cutoff = 0.001; % cutoff for vacuum
    %contsteps = 1000;
    
    %% GPU data transfers %%
    % n = gpuArray(n);
    % eps = gpuArray(eps);
    % D1x = gpuArray(D1x);
    % D1y = gpuArray(D1y);
    % laplace = gpuArray(laplace);
    % Vhat = gpuArray(Vhat);
    % uvec = gpuArray(uvec);
    % tangentvec = gpuArray(tangentvec);
    % ntol = gpuArray(ntol);
    % ngmrestol = gpuArray(ngmrestol);
    % nnitermax = gpuArray(nnitermax);
    % nminstep = gpuArray(nminstep);
    % cutoff = gpuArray(cutoff);
    % contsteps = gpuArray(contsteps);
    % errorsinarow = zeros(1,1,"gpuArray");
    % mu_list = zeros(1,contsteps,"gpuArray");
    % vac_list = zeros(1,contsteps,"gpuArray");
    
    % preconditioner for gmres
    pc= @(duvec) [reshape(ifft2(fft2(reshape(duvec(1:end-4),n,n))./(laplace-1),'symmetric'),n*n,1);duvec(end-3:end)]; % preconditioner

    % do first newton
    duvec=uvec;
    u0vec=uvec;
    [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x,D1y,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
    %if verbose
    %    disp(['first solution computed, mu=' num2str(uvec(end)) ', norm=' num2str(norm(uvec(1:end-4),'inf'))])
    %end
    uoldvec=uvec; % save first solution
    
    % do second newton
    uvec=uvec+ds*tangentvec; % take step
    u0vec=uvec;
    [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x,D1y,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep);
    %if verbose
    %    disp(['second solution computed, mu=' num2str(uvec(end)) ', norm=' num2str(norm(uvec(1:end-4),'inf'))])
    %end
    tangentvec=uvec-uoldvec;tangentvec=tangentvec/norm(tangentvec);
    
    
    h=figure(1);
    %w=waitbar(0,'forming vacuum...')
    
    % vac_pts = 0;
    % vac_pts_max = 12;
    
    %start secant continuation
    vactrue = 'f';
    fprintf(['computing iteration ' num2str(1,'%04u') ' of ' num2str(contsteps,'%04u') '... mu=' num2str(uvec(end),'%06.5f') ' umax=' num2str(norm(uvec(1:end-4),'inf'),'%06.5f') ' (vac: ' vactrue ')'])
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
            contcount = contcount-1;
            errorsinarow(1) = errorsinarow(1) + 1;
            if errorsinarow(1) == 10
                disp(['ERROR: continuous failure to converge'])
                break
            end    
        else
            errorsinarow(1) = 0;
            if verbose
                %disp([num2str(newtonflag.iter) ' newton iterations'])
                %disp([num2str(newtonflag.nresidual) ' residual after newton'])
		        %disp(['mu=' num2str(uvec(end)) ' umax=' num2str(norm(uvec(1:end-4),'inf'))])
            end
            tangentvec=uvec-uoldvec;tangentvec=tangentvec/norm(tangentvec);
            set(0,'CurrentFigure',h)
                u = reshape(uvec(1:end-4),n,n);
                surf(x,y,u,LineStyle="none");
                title(['$\mu=$' num2str(uvec(end)) ', $|u|_\infty$=' num2str(norm(uvec(1:end-4),'inf')) ' $ds=$' num2str(ds) ' a=' num2str(sum(uvec(1:end-4)))],'interpreter','latex')
                axis([-pi pi -pi pi 0.001 2.5])
                xlabel('x')
                ylabel('y')
                zlabel('u')
                drawnow
            mu_list(contcount) = uvec(end);
            % inf_list(contcount) = norm(uvec(1:end-4),'inf');
    
            % calculate size of vacuum
            save("test",'uvec');
            vacuum = uvec(1:end-4) < cutoff;
            vacuum_size = sum(vacuum)*4*pi*pi/(n*n);
            if vacuum_size > 0
		        disp(['vacuum reached after ' num2str(contcount) ' iterations'])    
		        break
                vactrue = 't';
            %     vac_pts = vac_pts+1;
            %     waitbar(vac_pts/vac_pts_max,w,'taking vacuum data...')
            end
            %vac_radius = sqrt(vacuum_size/(2*pi));
            %if verbose
                %disp(['vacuum size...' num2str(vacuum_size)])
            %end    
            vac_list(contcount) = vacuum_size;
            %l_list = [l_list;vac_radius];
            
            fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' num2str(contcount+1,'%04u') ' of ' num2str(contsteps,'%04u') '... mu=' num2str(uvec(end),'%06.5f') ' umax=' num2str(norm(uvec(1:end-4),'inf'),'%06.5f') ' (vac: ' vactrue ')'])
            if newtonflag.iter >5
                ds=ds/2; % decrease step size if newton not converging
            end
            if newtonflag.iter <3
                ds=min(1,ds*1.2); % increase stepsize if newton performing well
            end
        end
        % if max(u(1:end-3))>2.3
        %     display("reached breakpoint")
        %     break
        % end
        % if vac_pts == vac_pts_max
        %     break
        % end
    end
    mu_list_out = mu_list %gather(mu_list);
    vac_list_out = vac_list %gather(vac_list);
end

% figure(33)
% plot(mu_list,inf_list,'.-')
% xlabel(['$\mu$'],'Interpreter','latex')
% ylabel(['$|u|_\infty$'],'Interpreter','latex')
% 
% 
% figure(34)
% plot(mu_list,l_list,'.-')
% xlabel(['$\mu$'],'Interpreter','latex')
% ylabel('l')
% hold on
% mu_list_pred = linspace(1/(pi*pi),max(mu_list),512);
% prediction = nthroot(2*pi*pi*pi*(mu_list_pred-1/(pi*pi)),4);
% plot(mu_list_pred,prediction,'.-')
% axis([0.1013 max(mu_list)+0.0001 0 max([max(l_list) max(prediction)+0.05])])
