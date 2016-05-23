example5a
example5b

subplot('Position',[0.525 0.95-0.4*Dy/Dx_*2 0.4 0.4*Dy/Dx_*2])
myy = cumtrapz(stot,stog);
myy2 = cumtrapz(stot2,stog2);

plot(stot/10,myy)
hold on
plot(stot2/10,myy2)
[c, in1] = min(abs(stot-7000));
[c, in2] = min(abs(stot2-2000));
plot(stot(in1)/10,myy(in1),'.')
plot(stot2(in2)/10,myy2(in2),'.')
ylabel('Strain')
xlabel('Time (s)')
bp = '~/Documents/MATLAB/activnet/data/';
cd(bp);
load('extendrec_meas')
% 
% go=[];
% stoss=[];
% stogs=[];
% stops=[];
% stot=[];
% for ind=1:size(allt,1)
%     mu = 100*abs(allp(ind,3));
% %     if(allp(ind,3)<0)
% %         mu=100*mu;
% %     end
%     tr = 1./allp(ind,10);
%     tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2);
%     sscale=10*(allp(ind,7)*mu)^(1/2)/allp(ind,5);
%     nscale = (allp(ind,2)/allp(ind,5))^2*allp(ind,6);
%     ntscale = allp(ind,2)^2*allp(ind,6)/allp(ind,5)/abs(allp(ind,3));
%     ntrue = nscale;
%     if(tr/ntscale<100)
%         ntrue = nscale*sqrt(tr/ntscale)/10;
%     end
%     tstop = find(allt(ind,:)==0,2);
%     if(length(tstop)>1)
%         tstop = tstop(2)-1;
%     else
%         tstop = length(allt(ind,:));
%     end
%     t = allt(ind,1:tstop);
%     sl = allf(ind,1:tstop);
%     g = allg(ind,1:tstop);
%     
%     if(allp(ind,6)==10)
%          stoss = [stoss mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))/sscale];
%          stogs = [stogs mean(g)];
%          stops = [stops mean(sl(find(abs(t-2*tr)==min(abs(t-2*tr))):end))];
%          stot = [stot tr/tscale];
%                 plot(t,cumtrapz(t,g,2),'DisplayName',['\tau = ' num2str(1/allp(ind,10)) ',  \xi = ' num2str(allp(ind,6)) ',  \upsilon = ' num2str(allp(ind,7))])
% %             plot(t/tscale,sl/sscale,'DisplayName',['\tau = ' num2str(1/allp(ind,10))])
%      hold on
%     end
% end
% figure
% plot(stot,stoss,'.')
% figure
% plot(stops,stogs,'.')

