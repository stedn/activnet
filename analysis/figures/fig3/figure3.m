bp = '/Users/wmcfadden/ext_gap_llc/';
cd(bp);
load('allmeas')

sig = abs(allp(:,11)./allp(:,2).^2.*allp(:,5).^2./allp(:,6)*10000);


figure
alltau = [];
allG = [];
go=[];
for ind=1:size(allt,1)
    scale=allp(ind,2)^2/allp(ind,5)/abs(allp(ind,3))*allp(ind,6);
    subplot(2,1,1)
    if(sig(ind)==0.1)
        subplot(2,1,2)
    end
%     subplot(2,1,floor(log10(sig(ind))+3))
    t = allt(ind,:);
    sl2 = cumtrapz(t,allg(ind,:),2);
    sl = diff(log(sl2),1,2)./diff(log(t),1,2);
    cutoff = find(sl==min(sl));
    modelFun =  @(p,x) 1-p(1)*exp(-x/p(2));plot(abs(allp(:,3)),abs(allp(:,11))./mean(allg(:,10),2),'.')
    [coefEsts,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(t(cutoff:end-1), sl(cutoff:end), modelFun, [1 scale]);
    if(~isempty(t(cutoff+find(sl(cutoff:end)>0.5,1)))&&sl2(end)<0.1&&allp(ind,6)>0.05)
        alltau = [alltau; t(cutoff+find(sl(cutoff:end)>0.5,1)) scale];
        go=[go; 1];
        plot(t,sl2,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6)))])
        % plot(t,sl2/abs(allp(ind,11))*allp(ind,6)*allp(ind,2)^2/allp(ind,5)^2,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6)))])
        allG = [allG sl2(min(sl)==sl)];
        %         plot(t(1:end-1)/scale,sl,'.','DisplayName',[num2str(allp(ind,2)) ' ' num2str(abs(allp(ind,3))) ' ' num2str(allp(ind,5)) ' ' num2str(allp(ind,6)) ' ' num2str(abs(allp(ind,11)/allp(ind,6))) ' ' num2str(sl2(end))])
        hold on
    else
        go=[go; 0];
    end
    
end

go=boolean(go);

figure
plot(allp(go,2)./allp(go,5),abs(allp(go,11))./mean(allg(go,end/2:end),2)./allp(go,6),'.')

figure
plot(abs(allp(go,5)),abs(allp(go,11))./allG'./abs(allp(go,3)),'.')

figure
plot(allp(go,2).^2./allp(go,5)./abs(allp(go,3)).*allp(go,6)/4,alltau(:,1),'.')