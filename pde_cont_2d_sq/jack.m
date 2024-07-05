function dF = jack(duvec,uvec,u0vec,tangentvec,D1x,D1y,laplace,eps,Vhat,n)
    u = reshape(uvec(1:end-4),n,n); %recover nxn matrices, excl. appended params
    u0 = reshape(u0vec(1:end-4),n,n);
    du = reshape(duvec(1:end-4),n,n);

    uhat = fft2(u);
    duhat = fft2(du);

    extraparams = uvec(end-3:end); %extract parameters
    sx=extraparams(1);
    sy=extraparams(2);
    alpha = extraparams(3);
    mu = extraparams(end);

    dextraparams = duvec(end-3:end); %extract parameters
    dsx=dextraparams(1);
    dsy=dextraparams(2);
    dalpha = dextraparams(3);
    dmu = dextraparams(end);
    
    %Fx = D1x.*fft2(u.*ifft2(D1x.*uhat-mu*D1x.*Vhat.*uhat,'symmetric'))
    dFx = D1x.*fft2(du.*ifft2(D1x.*uhat-mu*D1x.*Vhat.*uhat,'symmetric'))+...           
            D1x.*fft2(u.*ifft2(D1x.*duhat-mu*D1x.*Vhat.*duhat,'symmetric'))+...
            D1x.*fft2(u.*ifft2(-dmu*D1x.*Vhat.*uhat,'symmetric')); % product rule
    dFy = D1y.*fft2(du.*ifft2(D1y.*uhat-mu*D1y.*Vhat.*uhat,'symmetric'))+...           
            D1y.*fft2(u.*ifft2(D1y.*duhat-mu*D1y.*Vhat.*duhat,'symmetric'))+...
            D1y.*fft2(u.*ifft2(-dmu*D1y.*Vhat.*uhat,'symmetric')); % product rule
    dF = dFx + dFy;
    

    dF0 = reshape(ifft2(eps*laplace.*duhat+sx*D1x.*duhat+dsx*D1x.*uhat+sy*D1y.*duhat+dsy*D1y.*uhat+dF,'symmetric')-alpha*du-dalpha*u,n*n,1);

    masscond = duhat(1,1); % mass condition
    
    duvectrimmed = duvec(1:end-4);
    u0vectrimmed = u0vec(1:end-4);
    %xphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed)
    xphasecond = duvectrimmed'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed); % phase conditions
    yphasecond = duvectrimmed'*(reshape(circshift(u0,1,2),n*n,1)-u0vectrimmed);

    secantcond = duvec'*tangentvec; % secant condition

    dF=real([dF0;masscond;xphasecond;yphasecond;secantcond]);
end