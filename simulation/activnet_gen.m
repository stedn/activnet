function p = activnet_gen(zet,L,mu,kap,lc,del,ups,phi,psi,r,sig,Dx,Dy,Df,Dw,ls,lf,tinc,tfin,seed)
    
%% this mess just ensures that any string input is converted to numbers
    if(ischar(zet)); zet = str2num(zet); end;
    if(ischar(L)); L = str2num(L); end;
    if(ischar(mu)); mu = str2num(mu); end;
    if(ischar(kap)); kap = str2num(kap); end;
    if(ischar(lc)); lc = str2num(lc); end;
    if(ischar(del)); del = str2num(del); end;
    if(ischar(ups)); ups = str2num(ups); end;
    if(ischar(phi)); phi = str2num(phi); end;
    if(ischar(psi)); psi = str2num(psi); end;
    if(ischar(r)); r = str2num(r); end ;
    if(ischar(sig)); sig = str2num(sig); end;
    if(ischar(Dx)); Dx = str2num(Dx); end;
    if(ischar(Dy)); Dy = str2num(Dy); end;
    if(ischar(Df)); Df = str2num(Df); end;
    if(ischar(Dw)); Dw = str2num(Df); end;
    if(ischar(ls)); ls = str2num(ls); end;
    if(ischar(lf)); lf = str2num(lf); end;
    if(ischar(tinc)); tinc = str2num(tinc); end;
    if(ischar(tfin)); tfin = str2num(tfin); end;
    if(ischar(seed)); seed = str2num(seed); end;
    
    ext = sig<0;
    
    rng(seed);
    
    %% use inputs to calculate number of filaments to add
    ncnt = ceil(L/ls)+1;
    N = floor(2*Dx*Dy/lc/L);
    
    nu=[];
    if(ups>0)
        nu = ups*double(rand(N,N)<phi).*(ones(N,N)-eye(N,N));
        nu = (nu+nu')/2;
    end
    
    %% initialize network
    
    p = zeros(N*ncnt,2);
    for i=1:N
        p((i-1)*ncnt+1,:) = [Dx*rand Dy*rand];
        thet = rand*2*pi;
        for j = 2:ncnt
            p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
        end
    end
    
    
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];

    
    fileID = 1;
    %% solve ode
    z0 = reshape(p,1,[]);

    if(tinc>0.05/2/r)
        tinc = 0.05/2/r;
    end
    tt = 0:tinc:tfin;
    
    fprintf(fileID,'%.3f',0);
    for i=1:length(z0)
        fprintf(fileID,' %.4f',z0(i));
    end
    fprintf(fileID,'\n');
    
    activnet(N,tt,z0,zet,L,mu,kap,del,nu,psi,sig,Dx,Dy,Df,Dw,ncnt,lf,ext,r,tinc,fileID);
    

end
