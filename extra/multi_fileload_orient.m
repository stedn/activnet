bp = '/Users/wmcfadden/xlrelax_all/';
cd(bp);
files = dir;
files = {files.name};
alltt = [];
alltest = [];
allstrain = [];
allstrun = [];
allvisc = [];
alldisp=[];
allenrg=[];
allangs=[];
details = [];
for f = files
    if(strfind(f{1},'_scr') )
            code = strsplit(f{1},'_');
            if(exist([code{1} '_out.txt'],'file'))
                code = code{1}
                fileload_orient
            end
    end
end

% bp = '/Users/wmcfadden/xlrelax_longtear/';
% cd(bp);
% files = dir;
% files = {files.name};
% for f = files
%     if(strfind(f{1},'_scr') )
%             code = strsplit(f{1},'_');
%             if(exist([code{1} '_out.txt'],'file'))
%                 code = code{1}
%                 fileload_orient
%             end
%     end
% end
% 
% bp = '/Users/wmcfadden/xlrelax_phase/';
% cd(bp);
% files = dir;
% files = {files.name};
% for f = files
%     if(strfind(f{1},'_scr') )
%             code = strsplit(f{1},'_');
%             if(exist([code{1} '_out.txt'],'file'))
%                 code = code{1}
%                 fileload_orient
%             end
%     end
% end
% 
figure
subplot(2,1,1)
xplot = allstrain(2:end,:);
yplot = abs(diff(allenrg)./diff(alltt));
plot((xplot(1:2:end-1,:)+xplot(2:2:end,:))/2,(yplot(1:2:end-1,:)+yplot(2:2:end,:))/2,'.')
% semilogx((allstrain(3:end,:)),alldisp(2:end,:))
subplot(2,1,2)
% semilogx((allstrain(3:end,:)),alltest(2:end,:))
yplot = alltest(2:end,:);
plot((xplot(1:2:end-1,:)+xplot(2:2:end,:))/2,(yplot(1:2:end-1,:)+yplot(2:2:end,:))/2,'.')


figure
subplot(2,1,1)
semilogx((allstrain(3:end,:)),alldisp(2:end,:).*repmat())
subplot(2,1,2)
semilogx((allstrain(3:end,:)),alltest(2:end,:))

figure
semilogx(exp(allstrain(2:end,:)),allenrg(2:end,:),'.')
hold on
semilogx(exp(allstrain(2:end-1,:)),alltest(2:end,:))

figure
plot3(allstrun(2:end-1,:),repmat(details(3,:)./details(6,:),size(alltest(2:end,:),1),1),alltest(2:end,:),'.')

