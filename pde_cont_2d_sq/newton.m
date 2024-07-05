function [uvec,newtonflag] = newton(uvec,u0vec,duvec,tangentvec,Vhat,D1x,D1y,laplace,eps,n,ntol,ngmrestol,pc,nnitermax,nminstep)
    %paramters for gmres
    restart = 200;
    ngmresitermax = 25;

    F = @(uvec) funk(uvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n);

    nrhs=F(uvec);           % compute residual
    nresidual=norm(nrhs);   % estimate residual
    %disp(['residual ' num2str(nresidual)])
    %pause

    newtonflag.error=0;
    iter=0;
    while nresidual>ntol
        dF = @(duvec) jack(duvec,uvec,u0vec,tangentvec,D1x,D1y,laplace,eps,Vhat,n);
        [nnincr,flag]=gmres(dF,nrhs,restart,ngmrestol,ngmresitermax,pc);
        if flag==1
            sprintf(['gmres did not converge, residual is ' num2str(nresidual) ' after ' num2str(iter) ' iterations'])
            newtonflag.error=1;
            break
        end

        uvec=uvec-nnincr; % Newton step
        iter=iter+1; % keep track of number of iterations

        % recompute residual
        nrhs=F(uvec);           % compute residual
        nresidual=norm(nrhs);   % estimate residual

        if iter>nnitermax
            sprintf(['Maximal number of newton iterations reached, giving up; residual is ' num2str(nresidual)])
            newtonflag.error=2;
            break
        end
        if  norm(nnincr)<nminstep
            sprintf(['Newton step is ineffective, giving up; residual is ' num2str(nresidual)])
            newtonflag.error=3;
            break
        end
    end

    newtonflag.iter=iter;
    newtonflag.nresidual=nresidual;
end
