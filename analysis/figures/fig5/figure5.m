example5a
example5b

subplot('Position',[0.525 0.94-0.4*Dy/Dx_*2 0.4 0.4*Dy/Dx_*2])
myy = cumtrapz(stot,stog);
myy2 = cumtrapz(stot2,stog2);

plot(stot/10,myy,'DisplayName','\tau_r = \infty')
hold on
ylabel('Strain')
xlabel('Time (s)')
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extendrec_meas')

indabl = find(allp(:,6)==1&allp(:,5)==0.5&allp(:,11)==-0.001&allp(:,12)==60);
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
    
    if(t(end)>1000)
                plot(t/10,cumtrapz(t,g,2),'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10)])
     hold on
    end
end
legend('Location','southeast')
[c, in1] = min(abs(stot-7000));
[c, in2] = min(abs(stot2-2000));
plot([0 1000],[0.4 0.4],':','Color',[0.25 0.25 0.25])
plot(stot(in1)/10,myy(in1),'.')
plot(stot2(in2)/10,myy2(in2),'.')
ylim([0 0.6])
annotation('textbox', [0.005 0.91 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.45 0.91 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


print('-depsc','-r0',['figure5a.eps']);

h2=figure;
example5c
example5d
example5e
example5f

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('actrec_meas')

subplot('Position',[0.525 topp-0.3 0.4 0.275])

indabl = find(allp(:,6)==10&allp(:,5)==0.3&allp(:,7)==1&allp(:,8)==0.25&...
    (allp(:,10)==0.0003|allp(:,10)==0.003|allp(:,10)==0.03|allp(:,10)==0.3|allp(:,10)==3));
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
                plot([0 t]/10,[0 sl],'DisplayName',['\tau_r = ' num2str(1/allp(ind,10)/10)])
                hold on
    end
end
ylabel('Stress (nN)')
xlabel('Time (s)')
legend('Location','northeast')
xlim([0 50])
annotation('textbox', [0.005 0.66 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.44 0.63 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


print('-depsc','-r0',['figure5b.eps']);

