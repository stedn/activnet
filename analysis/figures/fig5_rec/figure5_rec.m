example5a
example5b


subplot('Position',[0.575 0.94-0.2*Dy/Dx_*2 0.4 0.2*Dy/Dx_*2])


hold on
ylabel('Strain')
xlabel('Time (s)')
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extendrec_meas')

indabl = find(allp(:,6)==1&allp(:,5)==0.5&allp(:,11)==-0.0002&allp(:,12)==75);
[dum,srt] = sort(allp(indabl,10));
for ind=indabl(srt)'
    
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
                plot(t/10,cumtrapz(t,g,2),'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) ])
     hold on
    end
end
legend('Location','northeast')



subplot('Position',[0.075+0.0125 0.94-0.25*2 0.375 0.2])

for ind=indabl(srt)'
    tr = 1/allp(ind,10)/10;
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
    sscale=10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    nscale = (allp(ind,2)/allp(ind,5))^2*allp(ind,6);
    
    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
   
    t = allt(ind,1:tstop);
    g = allg(ind,floor(tstop/2):tstop);
    
    if(1)
         semilogx(tr/ntscale,abs(allp(ind,11))/mean(g),'.','DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) ])
         hold on
    end
end

ylabel('Effective Viscosity (nNs/\mum)')
xlabel('Normalized Recycling Time (\tau_r/\tau_c)')


load('extendrec_meas')

subplot('Position',[0.575+0.0125 0.94-0.25*2 0.375 0.2])
for ind=1:size(allt,1)
    
       
    mu = 100*abs(allp(ind,3));
    
    tr = 1./allp(ind,10);
    nscale = pi/4*(allp(ind,2)/allp(ind,5)-1)^2*allp(ind,6);
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
    
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    if(allp(ind,11)<0&&tr>=allp(ind,6))
                loglog(tr/ntscale,abs(allp(ind,11))/mean(g(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/nscale/2,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \sigma = ' num2str(allp(ind,11))])
     hold on
    end
end

load('extend_meas2')

for ind=1:size(allt,1)
    mu = abs(allp(ind,3));
    
    nscale = pi/4*(allp(ind,2)/allp(ind,5)-1)^2*allp(ind,6);
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
    
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    sl2 = cumtrapz(t,g,2);
    sl = diff(log(sl2),1,2)./diff(log(t),1,2);

    cutoff = find(sl==min(sl),1);
    term = cutoff+find(sl2(cutoff:end)>0.15,1);
    if(isempty(term))
        term = length(sl);
    end
    spt = cutoff+find(sl(cutoff:term-1)>0.75,1);
    eta=abs(allp(ind,11))./mean(allg(ind,spt:end),2);
    if(~isempty(spt))
                loglog(10^4.7+100*randn,eta/nscale,'o','Color',[0.25 0.25 0.25],'MarkerSize',6,'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
     hold on
    end
end
myx = logspace(-5,5,35);
loglog(myx,1./(1+1./myx.^(3/4)),'--')

loglog([0.00001 100000],[1 1],':','Color',[0.25 0.25 0.25])
ylim([0.001 5])
ylabel('Normalized Viscosity (\eta/\eta_c)')
xlabel('Normalized Recycling Time (\tau_r/\tau_c)')

annotation('textbox', [0.00 0.915 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.49 0.915 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.00 0.61 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.47 0.61 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


print('-depsc','-r0',['figure5a.eps']);



