%% figure4a
example4a
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('contract_meas')

ax1 = subplot('Position',[0.565,0.915-exw*Dy/Dx_-0.225,0.39,0.215]);
unis = unique(-allp(:,7)./allp(:,3));
cols = parula(length(unis)+1);
for muc = unique(-allp(:,3))'
    for ups = unique(allp(:,7))'
        subind = allp(:,3)==-muc&allp(:,7)==ups&allp(:,6)==10;
        coli = find(unis==ups/muc);
        x = -allp(subind,end);
        if(length(x)>2)
            w0 = repmat(allw(subind,1),1,size(allw(subind,:),2));
            y = max(-(allw(subind,:)-w0)./w0,[],2);
            [sortedX, sortIndex] = sort(x);
            plot(sortedX,y(sortIndex),'DisplayName',['\upsilon/\mu = ' num2str(ups/muc)],'Color',cols(coli,:))
            hold on
        end
    end
end
colormap(ax1,cols)
cb = colorbar('Location','south','Ticks',linspace(0,1,5),'TickLabels',{'0.3','1.0','3.3','10','33'});
title(cb,'\upsilon/\mu_c');
pos = cb.Position;
cb.Position = [pos(1)+pos(3)/2 pos(2) pos(3)/2 pos(4)/3];
xlabel('Stiffness Asymmetry (\mu_e/\mu_c)')
ylabel('Maximum Strain')

annotation('textbox', [0.01 0.89 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.01 0.68 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.68 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


print('-depsc','-r0',['figure4a.eps']);

%% figure4c

example4b


bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('astress_meas')

ax2 = subplot('Position',[0.565,0.915-exw*Dy/Dx_-0.225,0.39,0.215]);
for ind=1:size(allt,1)
    t = allt(ind,:);
    mx = max(allf(ind,:));
    if(mx>0.005&&allp(ind,end)>0)
        taui = find(allf(ind,:)>0.75*mx);
        stotau = t(taui(1));
        loglog((1/allp(ind,5)*allp(ind,6)/10/sqrt(allp(ind,7)*abs(100*allp(ind,3)))),stotau/10,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end
loglog([0.1,1000],[0.1,1000],'k:')
xlim([0.1 1000])
ylim([0.1 1000])
xlabel('$$\tau_a = $$ ($$\xi/l_c\sqrt{\mu_e\upsilon}$$)','interpreter','latex')
ylabel('Time of Max Stress')

annotation('textbox', [0.01 0.89 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.00 0.68 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.68 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0',['figure4b.eps']);





