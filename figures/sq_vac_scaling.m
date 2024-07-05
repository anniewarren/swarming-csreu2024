% plots comparison from continuation and theory for precoefficient and
% bifurcation onset
load('sq_vac_scaling.mat') % produced using pde_cont_2d_sq/gpu(cpu)_script
%close(23)
figure(23)
plot(x01n0128eps10000mulist,(x01n0128eps10000vaclist./2).^2,'Color','b','LineWidth',1)
hold on
plot(x02n0128eps01000mulist,(x02n0128eps01000vaclist./2).^2,'Color','m','LineWidth',1)
plot(x01n0512eps00500mulist,(x01n0512eps00500vaclist./2).^2,'Color','g','LineWidth',1)
%text(x01n0128eps10000mulist(end),0.75,'$\varepsilon$ = 0.1','Interpreter','latex')
%text(x02n0128eps01000mulist(end),0.75,'$\varepsilon$ = 0.01','Interpreter','latex')
%text(x01n0512eps00500mulist(end),0.7,'$\varepsilon$ = 0.005','Interpreter','latex')

% theory
mugrid = linspace(0.1,0.118,200);
plot(mugrid,2*pi^5*(mugrid-1/(pi^2)),'Color','k','LineWidth',1)
%text(0.10208,0.65,'thy','Interpreter','latex','HorizontalAlignment','right')

axis([0.1 0.118 0 0.8])
xlabel('$\mu$','Interpreter','latex')
ylabel('Area$^2$','Interpreter','latex')

legend({'$\varepsilon$ = 0.1','$\varepsilon$ = 0.01','$\varepsilon$ = 0.005','thy: $A^2 \approx 2\pi^5\mu$'},'Location','north','Interpreter','latex','Direction','reverse')
%legend('boxoff')

%saveas(23,'sq_vac_scaling')
%saveas(23,'sq_vac_scaling','jpg') % eps for importing to latex
