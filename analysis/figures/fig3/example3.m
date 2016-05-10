bp = '/Users/wmcfadden/extend_llc_ver/';
code = 'mtazskxf';%gcqbbcyr
cd(bp)

%% load param file and decipher params
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
zet=str2num(paree{2});L=str2num(paree{3});mu=-str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
xi=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});Dx=str2num(paree{13});Dy=str2num(paree{14});Df=str2num(paree{15});
Dw=str2num(paree{16});ls=str2num(paree{17});lf=str2num(paree{18});


%% load simulation data
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

%% store initial positions and initial measurements (all 0)
op = reshape(zt(1,:),[],2);
tl=0;
stof = 0;
stog = 0;
stoa = 0;
stot = 0;

   
%% setup timepoints and space points to measure
inds = 1:ceil(size(zt,1)/200):size(zt,1);
inds = inds(2:end);
bpos = linspace(0,Dx,51);
bpos = bpos(1:end-1)+bpos(2)/2;
ll = 4;
rl = 16;

%% for loop over timepoints to display

h1 = figure; 
clear mov mov2
h = figure('Position', [50, 100, 100+600*Dx/Dy, 600]);

lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
temp=flipud(winter(lst));
temp2 = bone(2*lst);
cc = (1-trp.^2).*temp2(lst+1:end,:)+(trp.^2).*temp(1:lst,:);
temp=hot(2*lst);
cc2 = (1-trp.^2).*temp2(lst+1:end,:)+(trp.^2).*temp(1:lst,:);
indi = 1;

for ind = inds
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    
    figure(h)
    clf
    netplot_str(p,L,lf,ls,Dx,Dy,cc,cc2,0.05);

    colormap([flipud(cc2);cc])
    colorbar('westoutside')
    drawnow
    mov(indi) = getframe(h);
    
    dp = (p-op);
    op = p;
    
    % remove data if has moved farther than realistically possible 
    % these events are due to crossing domain boundary or recycling
    jumpcut = 50*median(abs(dp(1:1:end,1)));
    subind = abs(dp(1:2:end-1,2))>jumpcut|abs(dp(2:2:end,2))>jumpcut|abs(dp(1:2:end-1,1))>jumpcut|abs(dp(2:2:end,1))>jumpcut;
    subind = [subind subind]';
    subind = subind(:);
    dp(subind,:)=[];
    p(subind,:)=[];
    
    % compute velocities
    v = dp/(t(ind)-tl);
    tl = t(ind);
    
    % compute filament strain
    [XY,sx,sy]=get_str(p,L,lf,ls,Dx,Dy);   

    % compute filament tension
    fx = mu*sx;
    fy = mu*sy;
    if(mu<0)
        fx = -fx.*(1+99*double(sx>0));
        fy = -fy.*(1+99*double(sy>0));
    end

    %bin tension data
    [bb,nb,sb]=bindata_line(XY,fx,bpos);
    [bc,nc,sc]=bindata_line(XY,abs(fx),bpos);
    
    
    subind = p(:,1)<=bpos(rl)&p(:,1)>=bpos(ll);
    
    % store data
    stot = [stot t(ind)];
    stog = [stog nanmean(v(subind,1)./(p(subind,1)-Dx*Dw))];
    stof = [stof nanmean(bb(ll:rl).*nb(ll:rl))/Dy];
    stoa = [stoa nanmean(bc(ll:rl).*nc(ll:rl))/Dy];
    
    % plot spatially resolved data 
    figure(h1)
    subplot(2,1,1)
    plot(p(~subind,1),v(~subind,1),'.')
    hold on
    plot(p(subind,1),v(subind,1),'.')
    hold off
    xlim([0,Dx/2])
    ylim([0,0.5*10^-3])
    subplot(2,1,2)
    plot(bpos,bb.*nb/Dy)
    hold on
    plot(bpos(ll:rl),bb(ll:rl).*nb(ll:rl)/Dy)
    hold off
    xlim([0,Dx/2])
    ylim([0,-1.1*sig])
    
    
    drawnow
    mov2(indi) = getframe(h1);
    
    indi=indi+1;
    
end

%% if there was any data to store we will now display it and save it
if(length(stof)>2)
    movie2avi(mov,[bp code '_mov_ex.avi']);
    close(h);
    movie2avi(mov2,[bp code '_data_ex.avi']);
    close(h1);
    
    h2 = figure;
    [ax,p1,p2] = plotyy(stot/10,stof,stot/10,cumtrapz(stot,stog));
    xlabel(ax(1),'Time (s)') % label x-axis
    ylabel(ax(1),'Stress (nN)') % label left y-axis
    ylabel(ax(2),'Strain') % label rigdat = [ht y-axis
    h_leg=annotation('textbox', [0.65 0.15 0.15 0.2],...
            'String',{['\xi = ' num2str(xi)],['L = ' num2str(L)],['l_c = ' num2str(lc)],...
            ['\mu = ' num2str(mu)],['\sigma = ' num2str(sig)]});
    print('-dpdf','-r0',[code '_ex.pdf']);
    close(h2)
    
    
end



