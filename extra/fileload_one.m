code = 'xzsqsmex';
A = importdata(['/Users/wmcfadden/xlrelax_better/' code '_out.txt']);
fid = fopen(['/Users/wmcfadden/xlrelax_better/' code '_scr.txt']);
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
t = A(:,1);
zt = A(:,2:end);
clear mov
figure('Position', [50, 100, 1200, 600]);
lst = size(zt,1);
trp = repmat((1:lst)'/lst,1,3);
cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
            
indi = 1;
for ind = 1:10:size(zt,1)
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

    netplot(p,L,lf,ls,D,cc(ind,:));
    
    mov(indi) = getframe;
    clf
    indi = indi +1;
end
movie2avi(mov,[code '_mov.avi']);

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