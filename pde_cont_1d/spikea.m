function F=spikea(u,uinit,tang,D1,D2,eps,V1h) 
    u0=u(1:end-3);
    s=u(end-2);
    alpha = u(end-1);
    mu = u(end);
    u0h = fft(u0);

    V_part_h = D1.*fft(u0.*ifft(D1.*u0h+mu*V1h.*u0h,'symmetric'));

    F0 = ifft(eps*D2.*u0h+s*D1.*u0h + V_part_h,'symmetric')-alpha*u0;
    F1 = u0h(1)-length(u0); % mass condition (first fourier coefficient is integral :)

    uinit0 = uinit(1:end-3);
    F2 = (u0-uinit0)'*(circshift(uinit0,1)-uinit0); % phase condition, roughly <u',uinit>
    F3 = (u-uinit)'*tang; % secant condition

    F=[F0;F1;F2;F3];
return
