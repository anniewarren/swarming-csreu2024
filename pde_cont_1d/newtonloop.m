function [u,newton_flag]=newtonloop(nresidual,ntol,du,u,uinit,tang,D1,D2,eps,V1h,nrhs,ngmrestol,pc,nnitermax,nminstep)


spike= @(u) spikea(u,uinit,tang,D1,D2,eps,V1h);

newton_flag.error=0;
nniter=0;
while (nresidual>ntol)
%
    dspike=@(du) dspikea(du,u,uinit,tang,D1,D2,eps,V1h);    % form Jacobian
    
    [nnincr,flag]=gmres(dspike,nrhs,200,ngmrestol,25,pc); 

                % gmres solve for increment
    if flag==1
        sprintf(['gmres did not converge, residual is ' num2str(nresidual) ' after ' num2str(nniter) ' iterations'])
        newton_flag.error=1;
        break
    end
    u=u-nnincr;                         % Newton step
    nniter=nniter+1;                % keep track of number of iterations
    %
    % recompute residual
    nrhs=spike(u);                     % compute residual
    nresidual=norm(nrhs) ;         % estimate residual
    %
    if nniter>nnitermax
        sprintf(['Maximal number of newton iterations reached, giving up; residual is ' num2str(nresidual)])
        
        newton_flag.error=2;
        break
    end
    %
    if  norm(nnincr)<nminstep
        sprintf(['Newton step is ineffective, giving up; residual is ' num2str(nresidual)])
        newton_flag.error=3;
        break
    end
%
end
newton_flag.nniter=nniter;
newton_flag.nresidual=nresidual;
end
