clear; close all;

%% 1D potential
mu = 1;
c = 1/8;
x = linspace(-pi,pi,1000);
potential = 1/(c*sqrt(2*pi))*exp(-x.^2/(2*c^2))-mu*cos(x);

figure(13)
plot(x,potential,LineWidth=2)
ax = gca;
ax.Visible = 'off';

%saveas(13,'genpotential1d','epsc')

%% 2D potential
a=1.4;
c=1/4;
figure(14)
n = 512;
x = linspace(-pi,pi,n)';
y = linspace(-pi,pi,n);
[X,Y] = meshgrid(x,y);
%xgrid = repmat(x,1,n);
%ygrid = repmat(y,n,1);
peaks = 1/(c*a)*exp(-(X.^2+Y.^2)/(2*c^2)) + 1/(c*a)*exp(-((X-pi).^2+(Y-pi).^2)/(2*c^2))...
    + 1/(c*a)*exp(-((X-pi).^2+(Y+pi).^2)/(2*c^2)) + 1/(c*a)*exp(-((X+pi).^2+(Y-pi).^2)/(2*c^2))...
    +1/(c*a)*exp(-((X+pi).^2+(Y+pi).^2)/(2*c^2));
potential3d = peaks-mu*(cos(X).*cos(Y));
surf(x,y,potential3d,EdgeColor="none")
axis([-pi pi -pi pi])

ax = gca;
ax.Visible = 'off';

%saveas(14,'genpotential2d','epsc')