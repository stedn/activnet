example4a

bp = '/Users/wmcfadden/activ_free_sweep_nonl/';
cd(bp);
load('allmeas')

ax1 = subplot('Position',[0.625,0.95-0.15,0.325,0.15]);
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
cb = colorbar('Location','south','Ticks',linspace(0,1,5),'TickLabels',{'0.3','1.0','3.3','10','33'});
title(cb,'\upsilon/\mu');
pos = cb.Position;
cb.Position = [pos(1)+pos(3)/2 pos(2) pos(3)/2 pos(4)/3];
xlabel('Stiffness Asymmetry (\mu_c/\mu_e)')
ylabel('Simulation Maximum Strain')

ax2 = subplot('Position',[0.625,0.95-0.35,0.325,0.15]);
for ind=1:size(allt,1)
    t = allt(ind,:);
    taui = find(allw(ind,:)>0.7*max(allw(ind,:)));
    stotau = t(taui(1));
    loglog(allp(ind,2)*allp(ind,6)./allp(ind,7),stotau/10,'.','DisplayName',[num2str(allp(ind,6))])
    hold on
end
loglog([1,1000],[1,1000],'k')
ylim([3 300])
xlim([3 300])
xlabel('Predicted \tau_m (L\xi/\upsilon)')
ylabel('Simulation Time of Max Strain')

% example4b

