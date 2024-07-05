% Generates main (background) plot
% load('bif_diag_inset_cont.mat')
% figure(51)
% plot(x01n0128eps01000mulist,x01n0128eps01000inflist,'.-')
%     %axis([0.007 0.015 1.1 2.2])
%     xlabel('$\mu$','Interpreter','latex')
%     ylabel('$|u|_\infty$','Interpreter','latex')
%     title('Bifurcation diagram with corresponding particle states')
% saveas(51,'bif_diag_main.fig')

% cleans inset plots to remove title/axes if needed
% fig = openfig('bif_diag_inset\n128_mu0.1040_iter006.fig','visible');
% ax = gca;
% title = get(ax,'title')
% title.Visible = "off";
% set(ax,'XColor', 'none','YColor','none')
% saveas(fig,'bif_diag_inset\n128_mu0.1040_iter006_clean.fig')

% fig = openfig('bif_diag_inset\n128_mu0.1024_iter006.fig','visible');
% ax = gca;
% title = get(ax,'title')
% title.Visible = "off";
% set(ax,'XColor', 'none','YColor','none')
% daspect([1 1])
% saveas(fig,'bif_diag_inset\n128_mu0.1024_iter006_clean.fig')

% open the first figure
fig1 = openfig('bif_diag_main.fig','visible');
hold on

%line(x,y) - make lines from diagram to little box (order left to right)
line([0.10401 0.10385], [2.1633 1.8],'Color','black')
line([0.10272 0.10325], [2.0001 1.77884],'Color','black')
line([0.10244; 0.103], [1.6845; 1.4],'Color','black')
%line([0.0466824; 0.048], [1.82246; 1.69],'Color','black','LineStyle','--')

% create the new axes
ax2 = axes(fig1,'Units','normalized','position',[0.6584 0.6322 0.15 0.15]);
set(ax2,'XColor', 'none','YColor','none')
% open the second figure
fig2 = openfig('bif_diag_inset\n128_mu0.1040_iter006_clean.fig','invisible');

% copy the stuff from the second figure's axes into the new axes
ch = findall(fig2.CurrentAxes);
ch(ch == fig2.CurrentAxes) = [];
copyobj(ch,ax2)
ax2.ZLim = [0.1 3];

% create the new axes
ax3 = axes(fig1,'Units','normalized','position',[0.4670 0.5775 0.15 0.15]);
set(ax3,'XColor', 'none','YColor','none')
% open the second figure
fig3 = openfig('bif_diag_inset\n128_mu0.1027_iter004_clean.fig','invisible');

% copy the stuff from the second figure's axes into the new axes
ch = findall(fig3.CurrentAxes);
ch(ch == fig3.CurrentAxes) = [];
copyobj(ch,ax3)
ax3.ZLim = [0.1 2];

ax4 = axes(fig1,'Units','normalized','position',[0.3542 0.3343 0.15 0.15]);
set(ax4,'XColor', 'none','YColor','none')
% open the second figure
fig3 = openfig('bif_diag_inset\n128_mu0.1024_iter003_clean.fig','invisible');

% copy the stuff from the second figure's axes into the new axes
ch = findall(fig3.CurrentAxes);
ch(ch == fig3.CurrentAxes) = [];
copyobj(ch,ax4)
ax4.ZLim = [0.1 2];

%saveas(fig1,'bif_diag_inset_cont.jpg')
            
% % create inset plots
% xpos = [0.6584; 0.4670; 0.2742];%positions for each little plot 
% ypos = [0.6322; 0.5775; 0.5143];