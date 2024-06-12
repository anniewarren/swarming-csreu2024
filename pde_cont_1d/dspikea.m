function dF=dspikea(du,u,uinit,tang,D1,D2,eps,V1h)
    u0=u(1:end-3);
    du0=du(1:end-3);
    s=u(end-2);
    alpha = u(end-1);
    kappa = u(end);
    ds=du(end-2);
    dalpha = du(end-1);
    dkappa = du(end);
    u0h = fft(u0);
    du0h = fft(du0);
    V_part_h = D1.*fft(du0.*ifft(D1.*u0h+kappa*V1h.*u0h,'symmetric'))+...           
            D1.*fft(u0.*ifft(D1.*du0h+kappa*V1h.*du0h,'symmetric'))+...
            D1.*fft(u0.*ifft(dkappa*V1h.*u0h,'symmetric')); % product rule
    dF0 = ifft(eps*D2.*du0h+s*D1.*du0h+ds*D1.*u0h+V_part_h,'symmetric')-alpha*du0-dalpha*u0;
    dF1 = du0h(1); % mass condition
    uinit0 = uinit(1:end-3);
    dF2 = du0'*(circshift(uinit0,1)-uinit0); % phase condition
    dF3 = du'*tang; % secant condition

    dF=real([dF0;dF1;dF2;dF3]);
return
