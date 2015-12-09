% bp = '/Users/wmcfadden/activ_free/';
% code = 'yzxfbrhp';
cd(bp)
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
% zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
% del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
% r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});
% Dx = 2*D;
% Dy = D;
zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
xi=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});Dx=str2num(paree{13});Dy=str2num(paree{14});Df=str2num(paree{15});
Dw=str2num(paree{16});ls=str2num(paree{17});lf=str2num(paree{18});

A = importdata([bp code '_out.txt']);
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
h = figure('Position', [50, 100, 100+600*Dx/Dy, 600]);
lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
temp=hot(2*lst);
cc2 = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*temp(1:lst,:);
edges = {linspace(0.5,1.5,50),linspace(-90,90,50)}  ;      
indi = 1;
dp = 0;
op = reshape(zt(1,:),[],2);
tl=0;
stof = [];
stoa = [];
stot = [];
inds = 2:floor(size(zt,1)/100):size(zt,1);%2:10:min(1000,size(zt,1));
for ind = inds
    stot = [stot t(ind)];
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    whitebg('black')
    set(gcf,'Color',[0 0 0])
    set(gcf,'InvertHardcopy','off')
    dp = (p-op);
    op = p;
    dp(dp(:,1)>Dx/4,1)=dp(dp(:,1)>Dx/4,1)-Dx;
    dp(dp(:,1)<-Dx/4,1)=dp(dp(:,1)<-Dx/4,1)+Dx;
    v = dp/(t(ind)-tl);
    tl = t(ind);
    bpos = linspace(0,Dx,11);
    bpos = bpos(1:end-1)+bpos(2)/2;
    [b,n,s]=bindata2(p(:,1),v(:,1),bpos);

    [XY,sx,sy]=get_str(p,L,lf,ls,Dx,Dy);   
    
    fx = mu*sx;
    if(mu<0)
        fx = -fx.*(1+99*double(sx>0));
    end
    
    fy = mu*sy;
    if(mu<0)
        fy = -fy.*(1+99*double(sy>0));
    end
    
    [bb,nb,s]=bindata_line(XY,fx,bpos);
    [bc,nc,s]=bindata_line(XY,abs(fx),bpos);
    [bb,nb,s]=bindata2((XY(:,1)+XY(:,3))/2,fx,bpos);
    [bc,nc,s]=bindata2((XY(:,1)+XY(:,3))/2,abs(fx),bpos);

    netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,0.1);


    axes('Position',[.75 .7 .12 .2])
    title('vel plot')
    box on
    plot(p(:,1),v(:,1),'y.');
    hold on;
    plot(bpos,b,'r','LineWidth',4);
    
    stof = [stof sum(fx)];
    stoa = [stoa sum(abs(fx))];
%     
    yex = 4*abs(sig)/xi*Dx/(L/lc)^2 ;
    plot(bpos,bb.*nb/max(bc.*nc)*yex/2,'b','LineWidth',3);
    plot(bpos,bc.*nc/max(bc.*nc)*yex/2,'g','LineWidth',3);
    plot(bpos,nc/max(nc)*yex/2,'p','LineWidth',3);
    ylim(yex*[-1 1])
    xlim([0,Dx])
    ylabel('velocity/stress/density')
    set(gca,'fontsize',14)
        h_leg=annotation('textbox', [0.75 0.2 0.12 0.45],'BackgroundColor',[1 1 1],...
            'String',{code,['xi = ' num2str(xi)],['L = ' num2str(L)],['lc = ' num2str(lc)],...
            ['mu = ' num2str(mu)],['sig = ' num2str(sig)],['phi = ' num2str(phi)],...
            ['r = ' num2str(r)],['stress = ' num2str(sum(abs(fx)))]});
        set(h_leg,'FontSize',16);
        set(h,'PaperPositionMode','auto')

    drawnow
    mov(indi) = getframe(h);
    clf
    indi = indi +1;
end
if(exist('mov'))
    movie2avi(mov,[bp code '_mov.avi']);
    h2 = figure;
    plot(stot,stof)
    hold on
    plot(stot,stoa)
    print('-dpng','-r0',[code '_fig.png']);
    close(h2)
end
close(h);