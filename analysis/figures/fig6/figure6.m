example6a
example6b

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('domain_meas')

subplot('Position',[exw+0.15 0.95-(0.8-exw)*Dy/Dx_*1.525 0.8-exw (0.8-exw)*Dy/Dx_*1.5])

indabl = find(allp(:,6)==10&allp(:,7)==0.1);
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
    
    if(allp(ind,6)==10&&allp(ind,7)==0.1&&t(end)>3000)
                plot(t/10,cumtrapz(t,g,2),'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10,4)])
     hold on
    end
end
legend('Location','northwest')
xlabel('Time (s)')
ylabel('Strain')


load('domain_meas')


indabl = find(allp(:,6)==10&allp(:,7)==0.1);
[dum,srt] = sort(allp(indabl,10));
subplot('Position',[0.1 0.9-0.3*2 0.3 0.22])

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
    
    if(allp(ind,6)==10&&allp(ind,7)==0.1&&t(end)>3000)
                semilogx(tr/tscale,mean(g(end-100:end)),'.')
     hold on
    end
end
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
ylabel('Strain Rate (1/s)','interpreter','latex')


subplot('Position',[exw+0.2 0.9-0.3*2 0.8-exw-0.07 0.25])

for ind=1:size(allp,1)
    mu = 100*abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
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
    
    if(t(end)>2*tscale)
         semilogx(tr/tscale,mean(g)/allp(ind,7)*allp(ind,6)*allp(ind,2),'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}L\xi/\upsilon$$)','interpreter','latex')

annotation('textbox', [0.005 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.005 0.9-0.4-0.01 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.42 0.9-0.4+0.01 0.05 0.05],'String','d)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


print('-depsc','-r0',['figure6a.eps']);






h2=figure;


load('domainllc_meas')


subplot('Position',[0.075 0.95-0.4 0.35 0.25])


for ind=1:size(allp,1)
    mu = 100*abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
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
    
    if(t(end)>tscale*10)
         plot(allp(ind,2),mean(g)*allp(ind,6),'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  L = ' num2str(allp(ind,2)) '  l_c = ' num2str(allp(ind,5)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Filament Length , L (\mum)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}\xi$$)','interpreter','latex')
ylim ([0 0.004])





subplot('Position',[0.575 0.95-0.4 0.35 0.25])
for ind=1:size(allp,1)
    mu = 100*abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
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
    
    if(t(end)>tscale*10)
         plot(allp(ind,5),mean(g)/allp(ind,7)*allp(ind,6)*allp(ind,2)^1.1,'.','Color',[0.25 0.25 0.25],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  L = ' num2str(allp(ind,2)) '  l_c = ' num2str(allp(ind,5)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Cross-link spacing, l_c (\mum)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}\xi L^{1.1} $$)','interpreter','latex')
ylim ([0 .11])

print('-depsc','-r0',['figure6b.eps']);
