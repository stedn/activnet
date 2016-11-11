%% figure2S
h2=figure;
orig_bp=pwd;
bp = '../../../data/';
cd(bp);
load('contract_meas')




subplot('Position',[0.075 0.95-0.2 0.35 0.2])
for ind=1:size(allt,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    
    t = allt(ind,:)/10;
    w = allw(ind,:);
    w = (w(1)-w)/w(1);
    taui = find(w>0.7*max(w));
    stotau = t(taui(1));
    loglog(L*xi/ups,stotau,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6))])
    hold on
end
loglog([0.2,200],[0.2,200],'k:')
xlim([0.2,200])
ylim([0.2,200])

xlabel('Predicted \tau_m = (L\xi/\upsilon)')
ylabel('Time of Max Strain')


load('astress_meas')

subplot('Position',[0.55 0.95-0.2 0.35 0.2])
for ind=1:size(allt,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    phi = allp(ind,8)
    
    mx = max(allf(ind,:));
    if(allp(ind,end)>0)
        plot(phi,mx*lc,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end

xlabel('Activity Fraction $$\phi$$','interpreter','latex')
ylabel('Normalized Max Stress ($$\sigma \cdot l_c$$)','interpreter','latex')


load('contract_meas')

ax1 = subplot('Position',[0.565,0.915-0.5-0.225,0.39,0.215]);
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




annotation('textbox', [0.02 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.42 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figureS2.eps']);
cd(orig_bp)

