bp = '/Users/wmcfadden/extend_rec_sweep/';
cd(bp);
load('allmeas')

sig = abs(allp(:,11)./allp(:,2).^2.*allp(:,5).^2./allp(:,6)*10000);


figure
alltau = [];
alleta = [];
allG = [];
go=[];
for ind=1:size(allt,1)
    scale=allp(ind,2)^2/allp(ind,5)/abs(allp(ind,3))*allp(ind,6);
    eta = (allp(ind,2)/allp(ind,5) - 1)^2*allp(ind,6);
    subplot(2,2,1)
    if(allp(ind,10)==0.01)
        subplot(2,2,2)
    end
    if(allp(ind,10)==0.1)
        subplot(2,2,3)
    end
    if(allp(ind,10)==1)
        subplot(2,2,4)
    end
%     subplot(2,1,floor(log10(sig(ind))+3))
    t = allt(ind,:);
    tstop = find(t==0,2);
    if(length(tstop)>1)
        tstop = tstop(2)-1;
    else
        tstop = length(t);
    end
    
    myt = t(1:tstop);
    myg = allg(ind,1:tstop);
    sl2 = cumtrapz(myt,myg,2);
    sl = diff(log(sl2),1,2)./diff(log(myt),1,2);
    cutoff = find(sl==min(sl));
    modelFun =  @(p,x) 1-p(1)*exp(-x/p(2));
%     [coefEsts,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(myt(cutoff:end-1), sl(cutoff:end), modelFun, [1 scale]);
    spt = cutoff+find(sl(cutoff:end)>0.8,1);
        
    if(allp(ind,11)<0&&~isempty(spt)&&allp(ind,6)>0.5)
        alltau = [alltau; t(cutoff+find(sl(cutoff:end)>0.5,1)) scale];
        go=[go; 1];
        spt = cutoff+find(sl(cutoff:end)>0.8,1);
        alleta = [alleta abs(allp(ind,11))./mean(allg(ind,spt:end),2)./eta];
%         plot(myt,sl2/eta,'.','DisplayName',[num2str(allp(ind,10)) ' ' num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(sig(ind))])
        % plot(t,sl2/abs(allp(ind,11))*allp(ind,6)*allp(ind,2)^2/allp(ind,5)^2,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6)))])
%         allG = [allG sl2(min(sl)==sl)];
                plot(myt(1:end-1)/scale,sl,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6))) ' ' num2str(sl2(end))])
        hold on
    else
        go=[go; 0];
    end
    
end
go = boolean(go)
figure
plot(1./(allp(go,10).*alltau(:,2)),alleta);