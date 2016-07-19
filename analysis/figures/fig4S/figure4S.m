%% figure4S
h2=figure;
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('contract_meas')




subplot('Position',[0.075 0.95-0.2 0.35 0.2])
for ind=1:size(allt,1)
    t = allt(ind,:);
    taui = find(allw(ind,:)>0.85*max(allw(ind,:)));
    stotau = t(taui(1));
    loglog(allp(ind,2)*allp(ind,6)/allp(ind,7),stotau/10,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6))])
    hold on
end
loglog([1,1000],[1,1000],'k:')
ylim([3 300])
xlim([3 300])
xlabel('Predicted \tau_m = (L\xi/\upsilon)')
ylabel('Time of Max Strain')


load('astress_meas')

subplot('Position',[0.55 0.95-0.2 0.35 0.2])
for ind=1:size(allt,1)
    t = allt(ind,:);
    mx = max(allf(ind,:));
    if(allp(ind,end)>0)
        plot(allp(ind,8),mx*allp(ind,5),'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end

xlabel('Activity Fraction $$\phi$$','interpreter','latex')
ylabel('Normalized Max Stress ($$\sigma \cdot l_c$$)','interpreter','latex')


% subplot('Position',[0.075 0.95-0.45 0.35 0.2])
% for ind=1:size(allt,1)
%     t = allt(ind,:);
%     mx = max(allf(ind,:));
%     if(mx>0.005&&allp(ind,end)>0)
%         plot(sqrt(allp(ind,7)*abs(100*allp(ind,3)))/allp(ind,5),mx,'k.','DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
%         hold on
%     end
% end
% 
% xlabel('Predicted $$\sigma_a$$ ($$\sqrt{\mu_e\upsilon}/l_c$$)','interpreter','latex')
% ylabel('Max Stress')

annotation('textbox', [0.02 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.42 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0',['figure4S.eps']);