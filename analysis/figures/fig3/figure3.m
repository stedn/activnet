example3

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extend_meas')

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
    term = cutoff+find(sl2(cutoff:end)>0.2,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.75,1);
    if(~isempty(spt))
        tau=myt(cutoff+find(sl(cutoff:term-1)>0.75,1));
        eta = abs(allp(ind,11))./mean(allg(ind,spt:end),2);
        axes(ax1);
        plot(allp(ind,2)./allp(ind,5)-1,eta/allp(ind,6),'.')
        hold on
        axes(ax2);
        if(~isempty(tau))
            loglog(allp(ind,2).^2./allp(ind,5)./abs(allp(ind,3)).*allp(ind,6)/10,tau/10,'.')
            hold on
        end
    end
    
end
axes(ax1);
ylim([0 10])
ylim([0 100])
xlabel('Number of crosslinks (L/l_c-1)')
ylabel('Normalized effective viscosity (\eta_{c}/\xi)')
plot(linspace(0,16,20),linspace(0,16,20).^2,'k:')
axes(ax2);
loglog([10 1000],[10 1000],'k:')
xlim([10 1000])
ylim([10 1000])
xlabel('Predicted \tau_c = (L^2\xi/l_c\mu)')
ylabel('Simulation Transition time ')
print('-depsc','-r0','figure3.eps');

