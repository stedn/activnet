example3
annotation('textbox', [0.015 0.93 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.440 0.93 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.440 0.70 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extend_meas2')
print('-depsc','-r0','figure3a.eps');
h2=figure;


ax1 = subplot('Position',[0.1,0.1,0.35,0.3]);
ax2 = subplot('Position',[0.55,0.1,0.35,0.3]);
for ind=1:size(allt,1)
    
    t = allt(ind,:);
    tstop = find(t==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(t);
    end
    
    myt = t(1:tstop);
    myg = allg(ind,1:tstop);
    sl2 = cumtrapz(myt,myg,2);
    sl = diff(log(sl2),1,2)./diff(log(myt),1,2);

    cutoff = find(sl==min(sl),1);
    term = cutoff+find(sl2(cutoff:end)>0.15,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.75,1);
    if(~isempty(spt))
        tau=myt(cutoff+find(sl(cutoff:term-1)>0.8,1));
        eta = abs(allp(ind,11))./mean(allg(ind,spt:end),2);
        axes(ax1);
        plot(allp(ind,2)./allp(ind,5)-1,eta/allp(ind,6),'.','Color',[0.25 0.25 0.25])
        hold on
        axes(ax2);
        if(~isempty(tau)&&tau<t(end))
            loglog(allp(ind,2).^2./allp(ind,5)./abs(allp(ind,3)).*allp(ind,6)/10,tau/10,'.','Color',[0.25 0.25 0.25])
            hold on
        end
    end
    
end
axes(ax1);
xlabel('Number of crosslinks per filament (L/l_c-1)')
ylabel('Normalized effective viscosity (\eta/\xi)')
plot(linspace(0,25,20),pi/4*linspace(0,25,20).^2,'k:')
axes(ax2);
loglog([10 1000],[10 1000],'k:')
xlim([10 1000])
ylim([10 1000])
xlabel('Predicted \tau_c (L^2\xi/l_c\mu_e)')
ylabel('Simulated Transition Time \tau_c (s)')

annotation('textbox', [0.005 0.37 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.46 0.37 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0','figure3b.eps');

