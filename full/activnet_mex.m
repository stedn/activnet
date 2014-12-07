%% setup initiation params
tfin = 100;
tinc = 0.1;

ncnt = 30;
lf = 0.01;
Df= 0.1;

zet = 0.1;
mu = 10;
kap = 0.01;

del = 10;
max_nu = 1;
fn = 1;

r = 0;

D = 4;
sig = 0;

L = 3;
lc = 0.2;
N = floor(4*D^2/lc/L);
    
nu=[];
if(max_nu>0)
    nu = max_nu*double(rand(N,N)<fn);
    nu = (nu+nu')/2;
end
%% initialize network
p = zeros(N*ncnt,2);
for i=1:N
    p((i-1)*ncnt+1,:) = D*[2*rand rand];
    thet = rand*2*pi;
    for j = 2:ncnt
        p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1)*[cos(thet) sin(thet)];
        thet = thet+pi/16*randn;
    end

end

p = [mod(p(:,1),2*D),mod(p(:,2),D)];


%% solve ode
z0 = reshape(p,[],1);
    
if(tinc<0.05*N*L/2/r)
    tinc = 0.05*N*L/2/r;
end
tt = 0:tinc:tfin;
zt = zeros(length(tt),length(z0));
zt(1,:)=z0';

tic

cfg = coder.config('mex');
argarg = {tfin,z0,zet,L,mu,kap,del,nu,sig,D,Df,ncnt,lf};
codegen -config cfg sp_activnet_mass.m -args argarg
codegen -config cfg activnet_ode.m -args argarg

options = odeset('Mass',@sp_activnet_mass_mex,'AbsTol',0.01,'RelTol',0.01);

figure('Position', [50, 100, 1200, 600]);
hold on
ind = 2;
while(ind<length(tt))
    [t,z] = ode23(@activnet_ode_mex,[tt(ind-1) tt(ind) tt(ind+1)],z0,options,zet,L,mu,kap,del,nu,sig,D,Df,ncnt,lf);
    zt(ind:ind+1,:) = z(2:3,:);
    
    %plot live
    p = reshape(zt(ind+1,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

    netplot(p,L,lf,ncnt,D);
    
    drawnow
    %end plot live
    
    z0 = z(3,:);
    if(r>0)
        i = randi(N,floor(r*2*tinc/L)+(rand<mod(r*2*tinc/L,1)),1);
        z0([(i-1)*ncnt+1; end/2+(i-1)*ncnt+1 ]) = D*[2*rand(size(i)); rand(size(i))];
        thet = rand*2*pi;
        for j = 2:ncnt
            z0([(i-1)*ncnt+j; end/2+(i-1)*ncnt+j ]) = z0([(i-1)*ncnt+j-1; end/2+(i-1)*ncnt+j-1 ])+L/(ncnt-1.0)*[cos(thet)*ones(size(i)); sin(thet)*ones(size(i))];
            thet = thet+pi/32*randn;
        end
    end
    
    ind = ind+2;
    
    toc
end


%% plot output
clear mov
figure('Position', [50, 100, 1200, 600]);
cc = copper(size(zt,1));
for ind = 1:size(zt,1)
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

    netplot(p,L,lf,ncnt,D);
    
    mov(ind) = getframe;
    clf
end

nbins = 30;
rs = linspace(0,1,nbins+1);
rs = rs(1:end-1);
sp = D*(rs+rs(2)/2);

figure;
xt = zt(:,1:end/2);
yt = zt(:,end/2+1:end);
for k=1:size(xt,1)
    sv = [];
    for r=rs
        subx = xt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
        suby = yt(:,max(xt(:,1:end/2))<(r+0.05)*D&min(xt(:,1:end/2))>r*D);
        sv = [sv nanmean(suby(k,:)-(suby(1,:)))];
    end
    ss(k)=mean(sv./sp);
end
plot(tt,ss)