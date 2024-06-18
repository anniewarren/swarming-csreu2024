function F = funk(uvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n)
    u = reshape(uvec(1:end-4),n,n); %recover nxn matrices, excl. appended params
    u0 = reshape(u0vec(1:end-4),n,n);

    uhat = fft2(u);

    extraparams = uvec(end-3:end); %extract parameters
    sx=extraparams(1);
    sy=extraparams(2);
    alpha = extraparams(3);
    mu = extraparams(end);
    
    %% F = div(u[grad(u)-mu(grad(V)*u)]) %%
    Fx = D1x.*fft2(u.*ifft2(D1x.*uhat-mu*D1x.*Vhat.*uhat,'symmetric'));
    Fy = D1y.*fft2(u.*ifft2(D1y.*uhat-mu*D1y.*Vhat.*uhat,'symmetric'));
    F=Fx+Fy;
    
    %%%%% need help here %%%%%
    F0 = reshape(ifft2(eps*laplace.*uhat+sx*D1x.*uhat + sy*D1y.*uhat + F,'symmetric')-alpha*u,n*n,1); %%% eps*D2*uhat, s*D1 to 2D??
    masscond = uhat(1,1)-n*n; % mass condition (first fourier coefficient is integral :)
    
    uvectrimmed = uvec(1:end-4);
    u0vectrimmed = u0vec(1:end-4);
    xphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed); % phase condition, roughly <u',uinit>
    yphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,2),n*n,1)-u0vectrimmed);
    secantcond = (uvec-u0vec)'*tangentvec; % secant condition

    F=[F0;masscond;xphasecond;yphasecond;secantcond];
end