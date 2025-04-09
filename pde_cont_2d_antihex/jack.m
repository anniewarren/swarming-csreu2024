function dF = jack(duvec,uvec,u0vec,tangentvec,D1y1,D1y2,laplace,eps,Vhat,n)
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

    D1x1 = D1y1;
    D1x2 = D1y1/sqrt(3) + 2*D1y2/sqrt(3);
    
    %Fx = D1x.*fft2(u.*ifft2(D1x.*uhat-mu*D1x.*Vhat.*uhat,'symmetric'))
    dFx1 = D1x1.*fft2(du.*ifft2(D1x1.*uhat-mu*D1x1.*Vhat.*uhat,'symmetric'))+...           
            D1x1.*fft2(u.*ifft2(D1x1.*duhat-mu*D1x1.*Vhat.*duhat,'symmetric'))+...
            D1x1.*fft2(u.*ifft2(-dmu*D1x1.*Vhat.*uhat,'symmetric')); % product rule
    dFx2 = D1x2.*fft2(du.*ifft2(D1x2.*uhat-mu*D1x2.*Vhat.*uhat,'symmetric'))+...           
            D1x2.*fft2(u.*ifft2(D1x2.*duhat-mu*D1x2.*Vhat.*duhat,'symmetric'))+...
            D1x2.*fft2(u.*ifft2(-dmu*D1x2.*Vhat.*uhat,'symmetric')); % product rule
    dF = dFx1 + dFx2;
    
    
    dF0 = reshape(ifft2(eps*laplace.*duhat+sx*D1x1.*duhat+dsx*D1x1.*uhat+sy*D1x2.*duhat+dsy*D1x2.*uhat+dF,'symmetric')-alpha*du-dalpha*u,n*n,1);

    masscond = duhat(1,1); % mass condition
    
    duvectrimmed = duvec(1:end-4);
    u0vectrimmed = u0vec(1:end-4);
    %xphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed)
    xphasecond = duvectrimmed'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed); % phase conditions
    yphasecond = duvectrimmed'*(reshape(circshift(u0,1,2),n*n,1)-u0vectrimmed);

    secantcond = duvec'*tangentvec; % secant condition

    dF=real([dF0;masscond;xphasecond;yphasecond;secantcond]);
end