

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
