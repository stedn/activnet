example6aS
example6bS

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);

annotation('textbox', [0.005 0.92 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.5 0.92 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0',['figure6S.eps']);
load('domain_meas')

h2=figure
for ind=1:size(allp,1)
    mu = 100*abs(allp(ind,3));
    alln{ind}
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
    sscale=1/10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    if(1)
                semilogx(tr/tscale,mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/sscale,'.','DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10,3)])
                hold on
    end
end