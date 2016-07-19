bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);

topp=0.95;
w=0.43;
h=0.2;
tr = logspace(0,6,50);







xspot=0.05;
subplot('Position',[xspot topp-h w h])
tc = 100;
ta = 1000;
nr = 1./(1+(tc./tr).^(3/4));
sr = 1./(tr/ta + ta./tr);

loglog(tr,sr,'DisplayName','Stress (\sigma_{ss})')
hold on
loglog(tr,nr,'DisplayName','Viscosity (\eta_c)')
ylim([0 1])
xlabel('\tau_r')
set(gca,'ytick',[],'xtick',[tc ta], 'xticklabel',{'\tau_c' '\tau_a'},'xminortick','off')
title('\tau_c<\tau_a')
legend('Location','northeast')

subplot('Position',[xspot topp-h*2-0.05 w h])
semilogx(tr,sr./nr,'Color',[0.25 0.25 0.25])

ylabel('Strain Rate (1/s)')
xlabel('\tau_r')
set(gca,'ytick',[],'xtick',[tc ta], 'xticklabel',{'\tau_c' '\tau_a'},'xminortick','off')



















xspot=xspot+0.48;
subplot('Position',[xspot topp-h w h])
tc = 10000;
ta = 1000;
nr = 1./(1+(tc./tr).^(3/4));
sr = 1./(tr/ta + ta./tr);

loglog(tr,sr)
hold on
loglog(tr,nr)
ylim([0 1])
xlabel('\tau_r')
set(gca,'ytick',[],'xtick',[ta tc], 'xticklabel',{'\tau_a' '\tau_c'},'xminortick','off')
title('\tau_a<\tau_c')

subplot('Position',[xspot topp-h*2-0.05 w h])
semilogx(tr,sr./nr,'Color',[0.25 0.25 0.25])

ylabel('Strain Rate (1/s)')
xlabel('\tau_r')
set(gca,'ytick',[],'xtick',[ta tc], 'xticklabel',{'\tau_a' '\tau_c'},'xminortick','off')



annotation('textbox', [0.01 0.9 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.9 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.01 0.9-h-0.05 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.9-h-0.05 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0',['figure_theor.eps']);
