cd(bp)
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
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
edges = {linspace(0.5,1.5,50),linspace(-90,90,50)}  ;      
indi = 1;
dp = 0;
op = reshape(zt(1,:),[],2);
tl=0;
stof = [0];
stog = [0];
stoa = [0];
stot = [0];
inds = 1:ceil(size(zt,1)/1000):size(zt,1);%2:10:min(1000,size(zt,1));
inds = inds(2:end);
h1 = figure;
    
for ind = inds
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    
    dp = (p-op);
    dp(dp(:,1)>Dx/4,1)=dp(dp(:,1)>Dx/4,1)-Dx;
    dp(dp(:,1)<-Dx/4,1)=dp(dp(:,1)<-Dx/4,1)+Dx;
    
    jumpcut = 0.1;
    
    subind = abs(dp(1:2:end-1,2))>jumpcut|abs(dp(2:2:end,2))>jumpcut|abs(dp(1:2:end-1,1))>jumpcut|abs(dp(2:2:end,1))>jumpcut;
    subind = [subind subind]';
    subind = subind(:);
    
    op = p;
    
    dp(subind,:)=[];
    p(subind,:)=[];
    
    stot = [stot t(ind)];
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

    subind = p(:,1)<Dx*abs(Df)&p(:,1)>Dx*abs(Df)/2;
    stog = [stog mean(v(subind,1)./(p(subind,1)-Dx*Dw))];
    subplot(2,1,1)
    plot(p(subind,1),v(subind,1),'.')
    subplot(2,1,2)
    plot(p(subind,1),v(subind,1)./(p(subind,1)-Dx*Dw),'.')
    drawnow
    stof = [stof sum(fx)];
    stoa = [stoa sum(abs(fx))];
    indi = indi +1;
    
    
end
close(h1)
if(length(stof)>2)
    h2 = figure;
    [ax,p1,p2] = plotyy(stot,stof/Dy/abs(Dx*Df),stot,cumtrapz(stot,stog),'loglog','loglog');
    xlabel(ax(1),'Time') % label x-axis
    ylabel(ax(1),'Stress/Length') % label left y-axis
    ylabel(ax(2),'Strain') % label rigdat = [ht y-axis
    h_leg=annotation('textbox', [0.15 0.7 0.2 0.2],'BackgroundColor',[1 1 1],...
            'String',{code,['xi = ' num2str(xi)],['L = ' num2str(L)],['lc = ' num2str(lc)],...
            ['mu = ' num2str(mu)],['sig = ' num2str(sig)],['r = ' num2str(r)]});
    set(h_leg,'FontSize',12);
    print('-dpng','-r0',[code '_fig.png']);
    close(h2)
    if(~isempty(allg))
        stog = [stog zeros(1,size(allg,2)-length(stog))];
        stof = [stof zeros(1,size(allg,2)-length(stof))];
        stot = [stot zeros(1,size(allg,2)-length(stot))];
        allg = [allg zeros(size(allg,1),length(stog)-size(allg,2))];
        allf = [allf zeros(size(allf,1),length(stof)-size(allf,2))];
        allt = [allt zeros(size(allt,1),length(stot)-size(allt,2))];
    end
    allg = [allg; stog];
    allf = [allf; stof];
    allt = [allt; stot];
    allp = [allp; zet L mu kap lc xi ups phi psi r sig Dx Dy Df Dw];
    alln = {alln{:} code};
end
% close(h);