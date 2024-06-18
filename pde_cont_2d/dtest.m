rng(17)
duvec=randn(n*n+4,1);duvec=duvec/norm(duvec,'inf');
%duvec=zeros(n*n+4,1);
%duvec(end)=1;

h=0.001;
res=funk(uvec+h*duvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n)-funk(uvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n)- jack(h*duvec,uvec,u0vec,tangentvec,D1x,D1y,laplace,eps,Vhat,n);   
max(abs(res))
figure(2)
plot(res)

h=h/10;
res=funk(uvec+h*duvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n)-funk(uvec,u0vec,tangentvec,Vhat,D1x,D1y,laplace,eps,n)- jack(h*duvec,uvec,u0vec,tangentvec,D1x,D1y,laplace,eps,Vhat,n);   
max(abs(res))
