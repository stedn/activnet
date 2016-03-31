% bp = '/Users/wmcfadden/extend_rec_sweep/';
% bp = '/Users/wmcfadden/extend_rec_sweep_b/';
bp = '/Users/wmcfadden/extend_llc_ver/';
cd(bp);
load('allmeas')

sig = abs(allp(:,11)./allp(:,2).^2.*allp(:,5).^2./allp(:,6)*10000);

figure
alltau = [];
alleta = [];
allG = [];
go=[];
for ind=1:size(allt,1)
    scale=allp(ind,2)^2/allp(ind,5)/abs(100*allp(ind,3))*allp(ind,6);
    eta = pi/4*(allp(ind,2)/allp(ind,5))^2*allp(ind,6);
%     subplot(2,2,1)
%     if(allp(ind,10)==0.01)
%         subplot(2,2,2)
%     end
%     if(allp(ind,10)==0.1)
%         subplot(2,2,3)
%     end
%     if(allp(ind,10)==1)
%         subplot(2,2,4)
%     end
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
    cutoff = 1;%find(sl<0.5,1);
%     modelFun =  @(p,x) 1-p(1)*exp(-x/p(2));
%     [coefEsts,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(myt(cutoff:end-1), sl(cutoff:end), modelFun, [1 scale]);
    spt = cutoff+find(sl(cutoff:end)>0.2,1);
    term = cutoff+find(sl2(cutoff:end)>2,1);
    if(isempty(term))
        term = length(sl);
    end
    if(allp(ind,11)<0&&~isempty(spt)&&term-spt>10&&allp(ind,6)>0.01&&allp(ind,2)==1)
        alltau = [alltau; t(spt) scale];
        go=[go; 1];
        alleta = [alleta abs(allp(ind,11))./mean(allg(ind,spt:term),2)./eta];
%         plot(myt,sl2/eta,'.','DisplayName',[num2str(allp(ind,10)) ' ' num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(sig(ind))])
        % plot(t,sl2/abs(allp(ind,11))*allp(ind,6)*allp(ind,2)^2/allp(ind,5)^2,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6)))])
%         allG = [allG sl2(min(sl)==sl)];
%         plot(myt(spt:term)/10,sl2(spt:term),'.','DisplayName',[alln{ind} ' ' num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6))) ' ' num2str(sl2(end)) ' ' num2str(abs(allp(ind,11))./mean(allg(ind,spt:end),2)./eta)])
%         plot(sl2(spt:term-1),sl(spt:term-1),'.','DisplayName',[alln{ind} ' ' num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6))) ' ' num2str(sl2(end)) ' ' num2str(abs(allp(ind,11))./mean(allg(ind,spt:end),2)./eta)])
        chnk = 10;
        xt = sl2(spt:term-1);
        xt = xt(1:chnk*floor(length(xt)/chnk));
        yt = abs(allp(ind,11))./diff(sl2(spt:term)).*diff(myt(spt:term))/eta;
        yt = yt(1:chnk*floor(length(yt)/chnk));
        x = mean(reshape(xt, chnk, []));
        y = mean(reshape(yt, chnk, []));
        plot(x,y,'DisplayName',[alln{ind} ' ' num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6))) ' ' num2str(sl2(end)) ' ' num2str(abs(allp(ind,11))./mean(allg(ind,spt:end),2)./eta)])
%         
        
        hold on
    else
        go=[go; 0];
    end
    
end
go = boolean(go);
figure
% plot(1./(allp(go,10).*alltau(:,2)),alleta,'.');
plot(9500*ones(size(alleta)),alleta,'.');