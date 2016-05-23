example6a
example6b
example6c

subplot('Position',[0.525 0.95-0.4*Dy/Dx_*2 0.4 0.4*Dy/Dx_*2])
myy = cumtrapz(stot,stog);
myy2 = cumtrapz(stot2,stog2);
myy3 = cumtrapz(stot3,stog3);

plot(stot/10,myy)
hold on
plot(stot2/10,myy2)
plot(stot2/10,myy3)
ylabel('Strain')
xlabel('Time (s)')

bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('domain_meas')

go=[];
stoss=[];
stogs=[];
stops=[];
stot=[];
for ind=1:size(allt,1)
    mu = 100*abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tr = 1./allp(ind,10);
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
    sscale=10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
    nscale = (allp(ind,2)/allp(ind,5))^2*allp(ind,6);
    ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
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
    sl = allf(ind,1:tstop);
    g = allg(ind,1:tstop);
    
    if(allp(ind,6)==1)
         stoss = [stoss mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/sscale];
         stogs = [stogs mean(g)];
         stops = [stops mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))];
         stot = [stot tr/tscale];
                plot(t,cumtrapz(t,g,2),'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
%             plot(t/tscale,sl/sscale,'DisplayName',['\tau = ' num2str(1/allp(ind,10))])
     hold on
    end
end
figure
plot(stot,stoss,'.')
figure
plot(stops,stogs,'.')

