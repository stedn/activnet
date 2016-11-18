orig_bp=pwd;
example4_1_thinning
example4_2_thinning


%% plot the strain curves for the recycling and the no recycling case
subplot('Position',[0.55 0.94-0.2*Dy/Dx_*2 0.4 0.2*Dy/Dx_*2])
plot(stot(2:end),ston(2:end))
hold on
plot(stot2(2:end),ston2(2:end))

hold on
ylabel('num filaments')
xlabel('Time (s)')






%% drop some nnotations
annotation('textbox', [0.00 0.915 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.915 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.00 0.61 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.47 0.61 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


cd('../../../figures')
print('-depsc','-r0',['figure4_thinning.eps']);
cd(orig_bp)