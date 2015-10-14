% bp = '/Users/wmcfadden/activ_begin/';
% code = 'qrzaabds';% txveifsg
cd(bp)
A = importdata([bp code '_out.txt']);
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});

% zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
% del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
% %         r=str2num(paree{10});sig=str2num(paree{11});D=str2num(paree{12});Df=str2num(paree{13});ls=str2num(paree{14});lf=str2num(paree{15});
% r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});
%         sig=str2num(paree{10});D=str2num(paree{11});Df=str2num(paree{12});ls=str2num(paree{13});lf=str2num(paree{14});
fclose(fid);
A = A.data;
if(size(A,1)==1)
    imp2 = importdata([bp code '_out.txt'],' ',9);
    if(isfield(imp2,'data'))
        A = [A;imp2.data];
    end
end
t = A(:,1);
zt = A(:,2:end);
clear mov
h = figure('Position', [50, 100, 1200, 600]);
lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
temp=hot(2*lst);
cc2 = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*temp(1:lst,:);
edges = {linspace(0.5,1.5,50),linspace(-90,90,50)}  ;      
indi = 1;
dp = 0;
op = reshape(zt(1,:),[],2);
tl=0
for ind = 2:ceil(size(zt,1)/200):size(zt,1)
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];
    whitebg('black')
    set(gcf,'Color',[0 0 0])
    set(gcf,'InvertHardcopy','off')
    dp = (p-op);
    op = p;
    dp(dp(:,1)>D/2,1)=dp(dp(:,1)>D/2,1)-2*D;
    dp(dp(:,1)<-D/2,1)=dp(dp(:,1)<-D/2,1)+2*D;
    v = dp/(t(ind)-tl);
    tl = t(ind);
    netplot_str(p,L,lf,ls,D,cc,cc2,0.1);
%     pa = p(1:2:end-1,:);
%     pb = p(2:2:end,:);
%     spots1 = pa(:,1)<D&pb(:,1)>D;
%     spots2 = pa(:,2)<D/2&pb(:,2)>D/2;
%     pa=[pa;pa(spots1,:)];
%     pb=[pb;[pb(spots1,1)-2*D,pb(spots1,2)]];
%     pa(spots1,1)=pa(spots1,1)+2*D;
%     pa=[pa;pa(spots2,:)];
%     pb=[pb;[pb(spots2,1),pb(spots2,2)-D]];
%     pa(spots2,2)=pa(spots2,2)+D;
%     po = atan2(pb(:,1)-pa(:,1),pb(:,2)-pa(:,2));
%     po(po>pi/2)=po(po>pi/2)-pi;
%     po(po<-pi/2)=po(po<-pi/2)+pi;
    
    axes('Position',[.75 .7 .12 .2])
    title('vel plot')
    box on
    bpos = linspace(0,2*D,11);
    bpos = bpos(1:end-1)+bpos(2)/2;
    [b,n,s]=bindata2(p(:,1),v(:,1),bpos);
    plot(p(:,1),v(:,1),'y.');
    hold on;
    plot(bpos,b,'r','LineWidth',4);
    bpos = linspace(0,2*D,11);
    bpos = bpos(1:end-1)+bpos(2)/2;
%     [b,nb,s]=bindata2(px,po,bpos);
%     [b,nt,s]=bindata2(px(abs(po)<pi/4),po(abs(po)<pi/4),bpos);
    
%     plot(bpos,2*nt./nb*max([abs(ups)/zet/del mean(abs(v))/2]),'g','LineWidth',3);
    ylim(2*max([abs(ups)/zet/del mean(abs(v))/2])*[-1 1])
    ylabel('velocity/order param')
    set(gca,'fontsize',14)
    h_leg=annotation('textbox', [0.75 0.2 0.12 0.45],'BackgroundColor',[1 1 1],...
        'String',{code,['zeta*del = ' num2str(zet*del)],['L = ' num2str(L)],['lc = ' num2str(lc)],['mu = ' num2str(mu)],['ups = ' num2str(ups)]});
    set(h_leg,'FontSize',16);
    set(h,'PaperPositionMode','auto')
    
    drawnow
    mov(indi) = getframe(h);
    clf
    indi = indi +1;
end
movie2avi(mov,[bp code '_mov.avi']);
close(h);
% nbins = 30;
% rs = linspace(0,1,nbins+1);
% rs = rs(1:end-1);
% sp = D*(rs+rs(2)/2);
% 
% figure;
% xt = zt(:,1:end/2);
% yt = zt(:,end/2+1:end);
% for k=1:size(xt,1)
%     sv = [];
%     for r=rs
%         subx = xt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
%         suby = yt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
%         sv = [sv nanmean(suby(k,:)-(suby(1,:)))];
%     end
%     ss(k)=-mean(sv./sp);
%     sr(k)=std(sp*ss(k)-sv);
% end
% xsi0 = 1;
% gam0 = 
% fito = fit(t,ss','(1-exp(-xsi*x))*gam','StartPoint', [xsi0, gam0]);