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


load('contract_meas')

ax1 = subplot('Position',[0.55 0.95-0.2 0.35 0.2])
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




load('astress_meas_correct')

ax2 = subplot('Position',[0.075 0.95-0.2*2-0.05 0.35 0.2])
for ind=1:size(allt,1)
    t = allt(ind,:)/10;

    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10

    mx = max(allf(ind,:));
    if(mx>0.005&&allp(ind,end)>0)
        taui = find(allf(ind,:)>0.9*mx);
        stotau = t(taui(1));
        tau_pred = L*xi/sqrt(ups*mu)
        loglog(tau_pred,stotau,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end

load('actrec_meas')

for ind=1:size(allt,1)
    t = allt(ind,:)/10;

    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10

    tscale=L*xi/(ups*mu)^(1/2)

    
    mx = max(allf(ind,:));
    if(mx>0.005&&allp(ind,end)>09&&tr/tscale>10)
        taui = find(allf(ind,:)>0.9*mx);
        stotau = t(taui(1));
        tau_pred = L*xi/sqrt(ups*mu)
        loglog(tau_pred,stotau,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end

loglog([0.1,1000],[0.1,1000],'k:')
xlim([0.1 1000])
ylim([0.1 1000])
xlabel('$$\tau_a = $$ ($$L\xi/\sqrt{\mu_e\upsilon}$$)','interpreter','latex')
ylabel('Time of Max Stress')




load('astress_meas_correct')

subplot('Position',[0.55 0.95-0.2*2-0.05 0.35 0.2])
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
    if(allp(ind,end)>0&&allp(ind,8)>0.1&&allp(ind,8)<0.9)
        loglog(0.1/2*ups^(2/3)*mu^(1/3)/lc,mx,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end


load('actrec_meas')

for ind=1:size(allt,1)
    r = allp(ind,10)*10
    tr = 1/r;
    L = allp(ind,2)
    lc = allp(ind,5)
    mu = 100*abs(allp(ind,3))
    xi = allp(ind,6)/10
    ups = allp(ind,7)
    phi = allp(ind,8)
    
    tscale=L*xi/(ups*mu)^(1/2)
    
    mx = max(allf(ind,:));
    if(allp(ind,end)>0&&allp(ind,8)>0.1&&allp(ind,8)<0.9&&tr/tscale>10)
        loglog(0.1/2*ups^(2/3)*mu^(1/3)/lc,mx,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
        hold on
    end
end

loglog([0.01/2,1],[0.01,1],'k:')
xlim([0.01/2 1])
ylim([0.01/2 1])
xlabel('Estimated Max Stress $$\upsilon^{2/3} \mu_e^{1/3}/l_c$$','interpreter','latex')
ylabel('Max Stress ($$\sigma$$)','interpreter','latex')


annotation('textbox', [0.02 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.45 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

cd('../figures')
print('-depsc','-r0',['figureS2.eps']);
cd(orig_bp)


% load('astress_meas')
% 
% subplot('Position',[0.565,0.915-0.5-0.225,0.39,0.215]);
% for ind=1:size(allt,1)
%     r = allp(ind,10)*10
%     tr = 1/r;
%     L = allp(ind,2)
%     lc = allp(ind,5)
%     mu = 100*abs(allp(ind,3))
%     xi = allp(ind,6)/10
%     ups = allp(ind,7)
%     phi = allp(ind,8)
%     
%     mx = max(allf(ind,:));
%     if(allp(ind,end)>0)
%         plot(phi,mx*lc,'.','Color',[0.25 0.25 0.25],'DisplayName',[num2str(allp(ind,6)) ' ' num2str(allp(ind,end))])
%         hold on
%     end
% end
% 
% xlabel('Activity Fraction $$\phi$$','interpreter','latex')
% ylabel('Normalized Max Stress ($$\sigma \cdot l_c$$)','interpreter','latex')
