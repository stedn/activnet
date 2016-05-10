bp = '/Users/wmcfadden/activ_free_sweep_nonl/';
cd(bp);
load('allmeas')


% plot(allt(end,2:end)/10,(allg(end,2:end)-allg(end,2))/allg(end,2))
% hold on
% plot(allt(end,2:end)/10,allc(end,2:end))
% plot(allt(end,2:end)/10,alle(end,2:end))

% figure
% stotau = [];
% stogam = [];
% for ind=1:size(allt,1)
% %     subplot(2,2,1)
% %     if(allp(ind,7)==0.1&&allp(ind,3)==-0.01)
% %         subplot(2,2,2)
% %     end
% %     if(allp(ind,7)~=0.1&&allp(ind,3)==-0.01)
% %         subplot(2,2,3)
% %     end
% %     if(allp(ind,7)~=0.1&&allp(ind,3)~=-0.01)
% %         subplot(2,2,4)
% %     end
%     t = allt(ind,2:end)-allt(ind,2);
%     st = -(allg(ind,2:end)-allg(ind,2))./allg(ind,2);
%     ct = find(st==max(st));
%     modelFun =  @(p,x) p(1)*(1-exp(-x/p(2)));
%     [coefEsts,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(t(1:ct), st(1:ct), modelFun, [0.5 100]);
%     ct = min(2*ct,length(st));
%     stogam = [stogam max(st)];
%     stotau = [stotau t(find(st==max(st)))];
%     nl = -allp(ind,end);
%     muc = -allp(ind,3);
%     mue = muc*nl;
%     ups = allp(ind,7);
%     xi = allp(ind,6);
%     muups_v = 1;
%     if(xi==10&&muc/ups>=muups_v*0.75&&muc/ups<=muups_v*1.25)
% %         subplot(3,1,1)
% %         semilogx(allt(ind,1:ct)*ups/xi,alle(ind,1:ct)*(mue/ups)^0.5,'DisplayName',[num2str(allp(ind,3)*allp(ind,end)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8)), ' ',  num2str(allp(ind,end))]);
% %         hold on;
% %         subplot(3,1,2)
% %         semilogx(allt(ind,1:ct)*sqrt(ups*muc)/xi,-allc(ind,1:ct)*sqrt(muc/ups),'DisplayName',[num2str(allp(ind,3)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8)), ' ',  num2str(allp(ind,end))]);
% %         hold on;
% %         subplot(3,1,3)
% %         semilogx(t(1:ct)*sqrt(ups*muc)/xi,st(1:ct)*sqrt(muc/ups),'DisplayName',[num2str(allp(ind,3)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8)), ' ',  num2str(allp(ind,end))]);
% %         hold on;
%         
%         
%         subplot(3,1,1)
%         semilogx(allt(ind,1:ct)*ups/xi,alle(ind,1:ct),'DisplayName',[num2str(allp(ind,3)*allp(ind,end)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8)), ' ',  num2str(allp(ind,end))]);
%         hold on;
%         subplot(3,1,2)
%         semilogx(allt(ind,1:ct)/xi*ups,-allc(ind,1:ct),'DisplayName',[num2str(allp(ind,3)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8)), ' ',  num2str(allp(ind,end))]);
%         hold on;
%         subplot(3,1,3)
%         semilogx(t(1:ct)/xi*ups,st(1:ct),'DisplayName',[num2str(allp(ind,3)), ' ',  num2str(allp(ind,6)), ' ',  num2str(allp(ind,7)), ' ',  num2str(allp(ind,8)), ' ',  num2str(allp(ind,end))]);
%         hold on;
%     end
% end
figure
for muc = unique(-allp(:,3))'
    for ups = unique(allp(:,7))'
%         for xi = unique(allp(:,6))'
            subind = allp(:,3)==-muc&allp(:,7)==ups&allp(:,6)==10;
            x = -allp(subind,end);
            y = stogam(subind);
            [sortedX, sortIndex] = sort(x);
            plot(sortedX,y(sortIndex),'DisplayName',['\mu = ' num2str(muc) ' \upsilon = ' num2str(ups)])
            hold on
%         end
    end
end
plot(allp(:,6)./allp(:,7),stotau/10,'.')
% subind = allp(:,6)==100&abs(allp(:,3))==0.001;
% plot(stogam(subind),stotau(subind),'.')
