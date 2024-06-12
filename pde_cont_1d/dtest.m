rng(1)
du=randn(n+3,1);du=du/norm(du,'inf');
% du=zeros(n+3,1);
du(end)=1;
h=1e-3;
res=spikea(u+h*du,uinit,tang,D1,D2,eps,V1h)-spikea(u,uinit,tang,D1,D2,eps,V1h) -dspikea(h*du,u,uinit,tang,D1,D2,eps,V1h);   
max(abs(res))
figure(2)
plot(res)

h=h/10;
res=spikea(u+h*du,uinit,tang,D1,D2,eps,V1h)-spikea(u,uinit,tang,D1,D2,eps,V1h) -dspikea(h*du,u,uinit,tang,D1,D2,eps,V1h);   
max(abs(res))
