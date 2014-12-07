%% setup initiation params
%compare to htiwdmfz
tfin = 10;
tinc = 0.1;

ls = 0.5;
lf = 0.01;
Df= 0.01;

zet = 0.005;
mu = 1;
kap = 0.000001;

del = 10000;
max_nu = 0;
phi = 0;
psi = 0;

r = 0;

D = 4;
sig = 0.2;

L = 2;
lc = 0.1;

ncnt = ceil(L/ls)+1;
N = floor(4*D^2/lc/L);

nu=[];
if(max_nu>0)
    nu = max_nu*double(rand(N,N)<phi);
    nu = (nu+nu')/2;
end
%% initialize network
% p = zeros(N*ncnt,2);
% for i=1:N
%     p((i-1)*ncnt+1,:) = D*[2*rand rand];
%     thet = rand*2*pi;
%     for j = 2:ncnt
%         p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1)*[cos(thet) sin(thet)];
%         thet = thet+pi/16*randn;
%     end
% 
% end

p = zeros(N*ncnt,2);
    for i=1:N
        p((i-1)*ncnt+1,:) = D*[2*rand rand];
        thet = rand*2*pi;
        for j = 2:ncnt
            p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
        end
    end

% p = zeros(2*ncnt,2);
% i=1;
% p((i-1)*ncnt+1,:) = [D/2-L/3 D/2];
% thet = 0;
% for j = 2:ncnt
%     p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1)*[cos(thet) sin(thet)];
% end
% i=2;
% p((i-1)*ncnt+1,:) = [D/2-L/4 D/2-L/4];
% thet = pi/2;
% for j = 2:ncnt
%     p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1)*[cos(thet) sin(thet)];
% end

p = [mod(p(:,1),2*D),mod(p(:,2),D)];


%% solve ode
z0 = reshape(p,[],1);
    
if(tinc>0.05/2/r)
    tinc = 0.05/2/r;
end
tt = 0:tinc:tfin;
zt = zeros(length(tt),length(z0));
zt(1,:)=z0';

tic

options = odeset('Mass',@sp_activnet_mass,'AbsTol',0.001,'RelTol',0.001);

figure('Position', [50, 100, 1200, 600]);
hold on
ind = 2;
while(ind<length(tt))
    [t,z] = ode15s(@activnet_ode_pull,[tt(ind-1) tt(ind) tt(ind+1)],z0,options,zet,L,mu,kap,del,nu,psi,sig,D,Df,ncnt,lf);
    zt(ind:ind+1,:) = z(2:3,:);
    
    %plot live
    p = reshape(zt(ind+1,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

    netplot(p,L,lf,ls,D,'b');
    
    drawnow
    %end plot livea
    
    ind = ind+2;
    z0 = z(3,:);
    if(r>0)
        i = randi(N,floor(r*2*tinc*N)+(rand<mod(r*2*tinc*N,1)),1);
        z0([(i-1)*ncnt+1; end/2+(i-1)*ncnt+1 ]) = D*[2*rand(size(i)); rand(size(i))];
        thet = rand*2*pi;
        for j = 2:ncnt
            z0([(i-1)*ncnt+j; end/2+(i-1)*ncnt+j ]) = z0([(i-1)*ncnt+j-1; end/2+(i-1)*ncnt+j-1 ])+L/(ncnt-1.0)*[cos(thet)*ones(size(i)); sin(thet)*ones(size(i))];
            thet = thet+pi/32*randn;
        end
    end
    toc
end


%% plot output
lst = 93;
clear mov
figure('Position', [50, 100, 1200, 600]);
trp = repmat((1:lst)'/lst,1,3);
cc = (1-trp).*(winter(lst)*0.75+0.25*spring(lst))+(trp).*copper(lst);
for ind = 1:7:lst
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

   netplot(p,L,lf,ls,D,cc(ind,:));
    drawnow
end

nbins = 30;
rs = linspace(0,1,nbins+1);
rs = rs(1:end-1);
sp = D*(rs+rs(2)/2);

figure;
xt = zt(:,1:end/2);
yt = zt(:,end/2+1:end);
for k=1:10
    sv = [];
    for r=rs
        subx = xt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
        suby = yt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
        sv = [sv nanmean(suby(k,:)-(suby(1,:)))];
    end
    ss(k)=mean(sv./sp);
end
plot(tt,ss)