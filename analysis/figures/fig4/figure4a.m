bp = '/Users/wmcfadden/activ_free_p/';
cd(bp);
load('allmeas')


% plot(allt(end,2:end)/10,(allg(end,2:end)-allg(end,2))/allg(end,2))
% hold on
% plot(allt(end,2:end)/10,allc(end,2:end))
% plot(allt(end,2:end)/10,alle(end,2:end))

figure
stotau = [];
stogam = [];
for ind=1:size(allt,1)
%     subplot(2,2,1)
%     if(allp(ind,7)==0.1&&allp(ind,3)==-0.01)
%         subplot(2,2,2)
%     end
%     if(allp(ind,7)~=0.1&&allp(ind,3)==-0.01)
%         subplot(2,2,3)
%     end
%     if(allp(ind,7)~=0.1&&allp(ind,3)~=-0.01)
%         subplot(2,2,4)
%     end
    t = allt(ind,2:end)-allt(ind,2);
    st = abs(allg(ind,2:end)-allg(ind,2))./allg(ind,2);
    ct = find(st==max(st));
    modelFun =  @(p,x) p(1)*(1-exp(-x/p(2)));
    [coefEsts,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(t(1:ct), st(1:ct), modelFun, [0.5 100]);
    stogam = [stogam coefEsts(1)];
    stotau = [stotau coefEsts(2)];
    if(allp(ind,3)==-0.01)
        plot(t(1:ct)/allp(ind,6).*sqrt(allp(ind,7)),st(1:ct)./sqrt(allp(ind,7)),'DisplayName',[num2str(allp(ind,3)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8))]);
        hold on;
    end
end

% subind = allp(:,6)==100&abs(allp(:,3))==0.001;
% plot(stogam(subind),stotau(subind),'.')
