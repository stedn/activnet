bp = '/Users/wmcfadden/extend_instatear/';
code = 'xnwcbuzb';%wthlwzyh
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
stof = [0];
stog = [0];
stoa = [0];
stot = [0];
inds = 1:floor(size(zt,1)/4):size(zt,1);%2:10:min(1000,size(zt,1));floor(size(zt,1)/4):size(zt,1)
inds = inds(2:end);
for ind = inds
    stot = [stot t(ind)];
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    
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
    
    whitebg('black')
    set(gcf,'Color',[0 0 0])
    set(gcf,'InvertHardcopy','off')
    netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,0.1);

    
    subind = p(:,1)<Dx*abs(Df)&p(:,1)>Dx*abs(Df)/2;
    stog = [stog mean(v(subind,1)./(p(subind,1)-Dx*Dw))];
    stof = [stof mean(bb(1:3).*nb(1:3))];
    stoa = [stoa sum(abs(fx))];
    
    colormap(cc);
    colorbar('westoutside')
    drawnow
    mov(indi) = getframe(h);
    clf
    indi = indi +1;
end
if(length(stof)>2)
    movie2avi(mov,[bp code '_mov_ex.avi']);
    h2 = figure;
    [ax,p1,p2] = plotyy(stot,stof/Dy,stot,cumtrapz(stot,stog));
    xlabel(ax(1),'Time') % label x-axis
    ylabel(ax(1),'Stress') % label left y-axis
    ylabel(ax(2),'Strain') % label right y-axis
    print('-dpng','-r0',[code '_fig_ex.png']);
    close(h2)
end
close(h);