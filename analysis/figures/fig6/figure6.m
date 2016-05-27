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
                plot(t/10,cumtrapz(t,g,2),'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10)])
     hold on
    end
end
legend('Location','northwest')
xlabel('Time (s)')
ylabel('Strain')
annotation('textbox', [0.005 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.39 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


print('-depsc','-r0',['figure6a.eps']);

h2=figure;
subplot('Position',[0.075 0.95-(0.8-exw)*Dy/Dx_*1.525-0.35 0.35 0.25])

for ind=1:size(allp,1)
    mu = 100*abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    sscale=10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    nscale = (allp(ind,2)/allp(ind,5))^2*allp(ind,6);
    ntrue = nscale;
    if(tr/ntscale<100)
        ntrue = nscale*sqrt(tr/ntscale)/10;
    end
    tstop = find(allt(ind,:)==0,1);
    if(length(tstop)>0)
        tstop = tstop(1)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    g = allg(ind,floor(tstop/2):tstop);
    
    if(t(end)>500)
         semilogx(tr/tscale,mean(g)/allp(ind,7)*allp(ind,6),'.','DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10) '  \xi = ' num2str(allp(ind,6)) '  \upsilon = ' num2str(allp(ind,7)) '  t_{last} = ' num2str(t(end))])
         hold on
    end
end
xlabel('Normalized Recycling Time (\tau_r/\tau_a)')
ylabel('Normalized Strain Rate ($$\dot{\gamma}\xi/\upsilon$$)','interpreter','latex')

print('-depsc','-r0',['figure6b.eps']);

