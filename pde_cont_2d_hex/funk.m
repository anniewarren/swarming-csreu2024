function F = funk(uvec,u0vec,tangentvec,Vhat,D1x1,D1x2,laplace,eps,n)
    u = reshape(uvec(1:end-4),n,n); %recover nxn matrices, excl. appended params
    u0 = reshape(u0vec(1:end-4),n,n);

    %%% plotting %%%
    if false
        L = pi;
        y1grid = linspace(-L,L,n+1)';y1grid=y1grid(1:end-1); %column vector in hex
        y2grid = linspace(-L,L,n+1);y2grid=y2grid(1:end-1); %row vector in hex
        [Y1,Y2] = meshgrid(y1grid,y2grid); % mesh in hex


        %%%%%%%%%%%%% change of basis components %%%%%%%%%%%%%
        % hex basis vectors
        e1=[1;0];
        e2=[-1/2;sqrt(3)/2];
        T = [e1,e2]; % square to hex change of basis matrix

        % calculating boundaries for plotting
        lowerLeftbound = T*[-L;-L];
        lowerRightbound = T*[L;-L];
        upperLeftbound = T*[-L;L];
        upperRightbound = T*[L;L];
        ylowlim = lowerRightbound(2);
        yuplim = upperRightbound(2);

        % axes for plotting hex
        X1 = T(1,1)*Y1 + T(1,2)*Y2;
        X2 = T(2,1)*Y1 + T(2,2)*Y2; % here T(2,1) = 0
        
        figure(13)
        surf(X1,X2,u,EdgeColor="none"); %3D surface plot
        axis([upperLeftbound(1) lowerRightbound(1) ylowlim yuplim 0.1 6])
        daspect([1 1 1])
        clim([0.1,6]);
        colorbar;
    end
    
    uhat = fft2(u);

    extraparams = uvec(end-3:end); %extract parameters
    sx=extraparams(1);
    sy=extraparams(2);
    alpha = extraparams(3);
    mu = extraparams(end);
    
    D1y1 = D1x1;
    D1y2 = D1x1/sqrt(3) + 2*D1x2/sqrt(3);

    %% F = div(u[grad(u)-mu(grad(V)*u)]) %%
    Fy1 = D1y1.*fft2(u.*ifft2(D1y1.*uhat-mu*D1y1.*Vhat.*uhat,'symmetric'));
    Fy2 = D1y2.*fft2(u.*ifft2(D1y2.*uhat-mu*D1y2.*Vhat.*uhat,'symmetric'));
    F=Fy1+Fy2;
   
    F0 = reshape(ifft2(eps*laplace.*uhat+sx*D1y1.*uhat + sy*D1y2.*uhat + F,'symmetric')-alpha*u,n*n,1); %%% eps*D2*uhat, s*D1 to 2D??
    masscond = uhat(1,1)-n*n; % mass condition (first fourier coefficient is integral :)
    
    uvectrimmed = uvec(1:end-4);
    u0vectrimmed = u0vec(1:end-4);
    xphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,1),n*n,1)-u0vectrimmed); % phase condition, roughly <u',uinit>
    yphasecond = (uvectrimmed-u0vectrimmed)'*(reshape(circshift(u0,1,2),n*n,1)-u0vectrimmed);
    secantcond = (uvec-u0vec)'*tangentvec; % secant condition

    F=[F0;masscond;xphasecond;yphasecond;secantcond];
end