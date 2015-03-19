bp = '/Users/wmcfadden/xlrelax_ext/';
code = 'zljqjtcz';
A = importdata([bp code '_out.txt']);
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
pare = strsplit(C{1}{10}, '>');
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
figure('Position', [50, 100, 1200, 600]);
lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
temp=hot(2*lst);
cc2 = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*temp(1:lst,:);
edges = {linspace(0.75,1.25,50),linspace(-90,90,50)}  ;      
indi = 1;
for ind = 1:100:size(zt,1)
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];
    whitebg('black')
    set(gcf,'Color',[0 0 0])
    set(gcf,'InvertHardcopy','off')

%     netplot_str(p,L,lf,ls,D,cc,cc2,0.25);
    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);
    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);
    subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)=subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)+2*D;
    subpL(subpL(:,2)<D/3&subpR(:,2)>2*D/3,2)=subpL(subpL(:,2)<D/3&subpR(:,2)>2*D/3,2)+D;
    subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)=subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)+2*D;
    subpR(subpR(:,2)<D/3&subpL(:,2)>2*D/3,2)=subpR(subpR(:,2)<D/3&subpL(:,2)>2*D/3,2)+D;
    subv = subpR-subpL;
    subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
    subpL = subpL;
    subpR = subpR;
    XY = [subpL subpR];
    str = sqrt((XY(:,3)-XY(:,1)).^2 + (XY(:,4)-XY(:,2)).^2);
    angs = atan((XY(:,4)-XY(:,2))./(XY(:,3)-XY(:,1)));
    subind = XY(:,1)<D*(1-Df)&XY(:,3)<D*(1-Df)&XY(:,1)>D*Df&XY(:,3)>D*Df;
    [N,C] = hist3([str(subind)/L,angs(subind)*180/pi],'Edges',edges);
    imagesc(C{2},C{1},N)
    mov(indi) = getframe;
    clf
    indi = indi +1;
end
movie2avi(mov,[code '_dist_mov.avi']);

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