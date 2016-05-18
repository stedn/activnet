% example3

bp = '/Users/wmcfadden/extend_llc_ver/';
cd(bp);
load('allmeas')

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
    spt = cutoff+find(sl(cutoff:end)>0,1);
    term = cutoff+find(sl2(cutoff:end)>0.2,1);
    if(isempty(term))
        term = length(sl);
    end
    if(~isempty(spt))
        tau=myt(spt);
        eta = abs(allp(ind,11))./mean(allg(ind,spt:end),2);
        axes(ax1);
        plot(allp(ind,2)./allp(ind,5),eta/allp(ind,6),'.')
        hold on
        axes(ax2);
        loglog(allp(ind,2).^2./allp(ind,5)./abs(allp(ind,3)).*allp(ind,6)/4,tau,'.')
        hold on
    end
    
end
axes(ax1);
ylim([0 100])
