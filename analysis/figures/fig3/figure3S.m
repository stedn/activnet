
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extend_meas2')
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
    spt = cutoff+find(sl(cutoff:term-1)>0.25,1);
    if(~isempty(spt)&&allp(ind,6)==1)
        tau=myt(cutoff+find(sl(cutoff:term-1)>0.8,1));
        G = abs(allp(ind,11))./sl2(spt);
        if(myt(spt)>500)
        axes(ax1);
        plot(myt(1:spt)/10,2*sl2(1:spt)/abs(allp(ind,11))*abs(allp(ind,3))/allp(ind,5))
        hold on
        end
        axes(ax2);
        plot(2*abs(allp(ind,3))/allp(ind,5),G,'.','Color',[0.25 0.25 0.25])
        hold on
        
    end
    
end
axes(ax2);
xlabel('Estimated elastic modulus (2\mu/l_c)')
ylabel('Measured elastic modulus, G_0 (nN/\mum)')
axes(ax1);
ylabel('Normalized Strain (\gamma/\sigma\cdot2\mu/l_c)')
xlabel('Time (s)')

annotation('textbox', [0.005 0.37 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.46 0.37 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0','figure3S.eps');

