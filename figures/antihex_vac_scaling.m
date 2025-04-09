load('antihex_vac_scaling.mat') % produced by pde_cont_2d_hex/script.m with desired parameters 
figure(24)
plot(x01n0128eps10000mulist,(x01n0128eps10000vaclist./4).^2,'.-','LineWidth',2)
hold on
plot(x02n0128eps01000mulist,(x02n0128eps01000vaclist./4).^2,'.-','LineWidth',2)
plot(x03n0256eps00500mulist,(x03n0256eps00500vaclist./4).^2,'.-','LineWidth',2)

%text(x01n0128eps10000mulist(end),0.75,'$\varepsilon$ = 0.1','Interpreter','latex')
%text(x02n0128eps01000mulist(end),0.75,'$\varepsilon$ = 0.01','Interpreter','latex')
%text(x03n0256eps00500mulist(end),0.7,'$\varepsilon$ = 0.005','Interpreter','latex')

mugrid = linspace(1/(sqrt(3)*pi*pi)-0.0001,0.7,512);
plot(mugrid,6*pi^5*(mugrid-1/(sqrt(3)*pi*pi)),'-','Color','k','LineWidth',2)
%text(1/(sqrt(3)*pi*pi),0.65,'thy','Interpreter','latex','HorizontalAlignment','right')

axis([1/(sqrt(3)*pi*pi)-0.0005 0.07 0 0.8])
xlabel('$\mu$','Interpreter','latex')
ylabel('Area$^2$','Interpreter','latex')

legend({'$\varepsilon$ = 0.1','$\varepsilon$ = 0.01','$\varepsilon$ = 0.005','thy: $A^2 \approx 6\pi^5\mu$'},'Location','north','Interpreter','latex','Direction','reverse')
legend('boxoff')

saveas(24,'antihex_vac_scaling')
saveas(24,'antihex_vac_scaling','pdf')
saveas(24,'antihex_vac_scaling','epsc')