function dF = jack(duvec,uvec,u0vec,tangentvec,D1x1,D1x2,laplace,eps,Vhat,n)
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

    D1y1 = D1x1;
    D1y2 = D1x1/sqrt(3) + 2*D1x2/sqrt(3);
    
    %Fx = D1x.*fft2(u.*ifft2(D1x.*uhat-mu*D1x.*Vhat.*uhat,'symmetric'))
    dFy1 = D1y1.*fft2(du.*ifft2(D1y1.*uhat-mu*D1y1.*Vhat.*uhat,'symmetric'))+...           
            D1y1.*fft2(u.*ifft2(D1y1.*duhat-mu*D1y1.*Vhat.*duhat,'symmetric'))+...
            D1y1.*fft2(u.*ifft2(-dmu*D1y1.*Vhat.*uhat,'symmetric')); % product rule
    dFy2 = D1y2.*fft2(du.*ifft2(D1y2.*uhat-mu*D1y2.*Vhat.*uhat,'symmetric'))+...           
            D1y2.*fft2(u.*ifft2(D1y2.*duhat-mu*D1y2.*Vhat.*duhat,'symmetric'))+...
            D1y2.*fft2(u.*ifft2(-dmu*D1y2.*Vhat.*uhat,'symmetric')); % product rule
    dF = dFy1 + dFy2;
    
    
    dF0 = reshape(ifft2(eps*laplace.*duhat+sx*D1y1.*duhat+dsx*D1y1.*uhat+sy*D1y2.*duhat+dsy*D1y2.*uhat+dF,'symmetric')-alpha*du-dalpha*u,n*n,1);

    masscond = duhat(1,1); % mass condition
    
    duvectrimmed = duvec(1:end-4);
    u0vectrimmed = u0vec(1:end-4);
    %xphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed)
    xphasecond = duvectrimmed'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed); % phase conditions
    yphasecond = duvectrimmed'*(reshape(circshift(u0,1,2),n*n,1)-u0vectrimmed);

    secantcond = duvec'*tangentvec; % secant condition

    dF=real([dF0;masscond;xphasecond;yphasecond;secantcond]);
end