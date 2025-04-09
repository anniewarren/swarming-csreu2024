function F = funk(uvec,u0vec,tangentvec,Vhat,D1y1,D1y2,laplace,eps,n)
    u = reshape(uvec(1:end-4),n,n); %recover nxn matrices, excl. appended params
    u0 = reshape(u0vec(1:end-4),n,n);
    
    uhat = fft2(u);

    extraparams = uvec(end-3:end); %extract parameters
    sx=extraparams(1);
    sy=extraparams(2);
    alpha = extraparams(3);
    mu = extraparams(end);
    
    D1x1 = D1y1;
    D1x2 = D1y1/sqrt(3) + 2*D1y2/sqrt(3);

    %% F = div(u[grad(u)-mu(grad(V)*u)]) %%
    Fx1 = D1x1.*fft2(u.*ifft2(D1x1.*uhat-mu*D1x1.*Vhat.*uhat,'symmetric'));
    Fx2 = D1x2.*fft2(u.*ifft2(D1x2.*uhat-mu*D1x2.*Vhat.*uhat,'symmetric'));
    F=Fx1+Fx2;
   
    F0 = reshape(ifft2(eps*laplace.*uhat+sx*D1x1.*uhat + sy*D1x2.*uhat + F,'symmetric')-alpha*u,n*n,1); %%% eps*D2*uhat, s*D1 to 2D??
    masscond = uhat(1,1)-n*n; % mass condition (first fourier coefficient is integral :)
    
    uvectrimmed = uvec(1:end-4);
    u0vectrimmed = u0vec(1:end-4);
    xphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed); % phase condition, roughly <u',uinit>
    yphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,2),n*n,1)-u0vectrimmed);
    secantcond = (uvec-u0vec)'*tangentvec; % secant condition

    F=[F0;masscond;xphasecond;yphasecond;secantcond];
end