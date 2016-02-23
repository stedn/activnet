bp = '/Users/wmcfadden/activ_rec_examples/';
cd(bp);
load('allmeas')



% figure
alltau = [];
alltaux=[];
allsig=[];
allsigx=[];
go=[];
for ind=1:size(allt,1)
    mu = abs(allp(ind,3));
%     if(allp(ind,3)<0)
%         mu=100*mu;
%     end
    tscale=allp(ind,6)/(allp(ind,7)*mu)^(1/2)*allp(ind,2)/allp(ind,5);
    sscale=(allp(ind,7)*mu)^(1/2)/2;
    tstop = find(allt(ind,:)==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(allt(ind,:));
    end
    t = allt(ind,1:tstop);
    sl = allf(ind,1:tstop);
    
    if(allp(ind,8)~=1&&allp(ind,3)<0)
        allsig = [allsig; max(sl)];
        alltau = [alltau; t(find(sl==max(sl),1))];
        allsigx= [allsigx; sscale];
        alltaux= [alltaux; tscale];
        
%         subplot(2,1,1)
%         if(allp(ind,5)==0.3)
%             subplot(2,1,2)
%         end
    %     modelFun =  @(p,x) 1-p(1)*exp(-x/p(2));
    %     [coefEsts,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(t(cutoff:end-1), sl(cutoff:end), modelFun, [1 scale]);
    %     if(allp(ind,7)==0.1)    
            plot(t,sl,'DisplayName',[alln{ind} num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,7))) ' ' num2str(abs(allp(ind,8)))])
    %     end
        %     plot(t(1:end-1)/scale,sl,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6)))])
        hold on
    end
end
% figure
% plot(allsigx,allsig,'.');
% figure
% plot(alltaux,alltau,'.');
% plot(sc2./sc3,mean(alla(subind,end-10:end),2),'.')
