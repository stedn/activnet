h2=figure;
example5c
example5d
example5e
example5f

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('actrec_meas')

subplot('Position',[0.525 topp-0.3 0.4 0.275])

indabl = find(allp(:,6)==10&allp(:,5)==0.3&allp(:,7)==0.1&allp(:,8)==0.25&...
    (allp(:,10)==0.0001|allp(:,10)==0.001|allp(:,10)==0.01|allp(:,10)==0.1));
[dum,srt] = sort(allp(indabl,10));
for ind=indabl(srt)'
    alln{ind}
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
                plot([0 t]/10,[0 sl],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10,4)])
                hold on
    end
end
ylabel('Stress (nN)')
xlabel('Time (s)')
legend('Location','northeast')
xlim([0 250])

indabl = find(allp(:,6)==10&allp(:,5)==0.3&allp(:,7)==0.1&allp(:,8)==0.25);
    %(allp(:,10)==0.0003|allp(:,10)==0.003|allp(:,10)==0.03|allp(:,10)==0.3|allp(:,10)==3));
[dum,srt] = sort(allp(indabl,10));
subplot('Position',[0.1 topp-0.3*2 0.3 0.2])
for ind=indabl(srt)'
    mu = 100*abs(allp(ind,3));
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
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
                semilogx(tr/tscale,mean(sl(end-10:end)),'.','DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10,3)])
                hold on
    end
end
ylabel('Steady State Stress (nN)')
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
xlim([0.0001 100])
set(gca,'XTick',[0.01 1 100],'XTickLabel',[0.01 1 100])

load('actrec_meas')


subplot('Position',[0.525+0.025 topp-0.3*2 0.4-0.05 0.23])
for ind=1:size(allt,1)
    mu = 100*abs(allp(ind,3));
    
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    tscale2=allp(ind,2)*allp(ind,6)/allp(ind,7);
    sscale=1/10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    if(1)
                loglog(tr/tscale,mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/sscale,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
     hold on
    end
end
myx = logspace(-4,4,35);
loglog(myx,1./(1./myx+myx),'--')
ylabel('Normalized Steady State Stress (\sigma/\sigma_a)')
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')

annotation('textbox', [0.005 0.66 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.66 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.001 0.27 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.3 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])

print('-depsc','-r0',['figure5b.eps']);

