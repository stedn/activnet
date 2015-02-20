%% setup initiation params  0.0005 9 1 0 0.25 100000 0 0 0 0 1.296 54 0.05 9 0.05 0.1 1000 1
function activnet_gen(zet,L,mu,kap,lc,del,ups,phi,psi,r,sig,D,Df,ls,lf,tinc,tfin,seed)
    if(ischar(zet))
        zet = str2num(zet);
    end
    if(ischar(L))
        L = str2num(L);
    end
    if(ischar(mu))
        mu = str2num(mu);
    end
    if(ischar(kap))
        kap = str2num(kap);
    end
    if(ischar(lc))
        lc = str2num(lc);
    end
    if(ischar(del))
        del = str2num(del);
    end
    if(ischar(ups))
        ups = str2num(ups);
    end
    if(ischar(phi))
        phi = str2num(phi);
    end
    if(ischar(psi))
        psi = str2num(psi);
    end
    if(ischar(r))
        r = str2num(r);
    end
    if(ischar(sig))
        sig = str2num(sig);
    end
    if(ischar(D))
        D = str2num(D);
    end
    if(ischar(Df))
        Df = str2num(Df);
    end
    if(ischar(ls))
        ls = str2num(ls);
    end
    if(ischar(lf))
        lf = str2num(lf);
    end
    if(ischar(tinc))
        tinc = str2num(tinc);
    end
    if(ischar(tfin))
        tfin = str2num(tfin);
    end
    if(ischar(seed))
        seed = str2num(seed);
    end
    ext = 0;
    if(sig<0)
        ext = 1;
    end
    
    rng(seed);
    
    %% convert inputs
    ncnt = ceil(L/ls)+1;
    N = floor(4*D^2/lc/L);
    
    nu=[];
    if(ups>0)
        nu = ups*double(rand(N,N)<phi).*(ones(N,N)-eye(N,N));
        nu = (nu+nu')/2;
    end
    
    %% initialize network
    p = zeros(N*ncnt,2);
    for i=1:N
        p((i-1)*ncnt+1,:) = D*[2*rand rand];
        thet = rand*2*pi;
        for j = 2:ncnt
            p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1.0)*[cos(thet) sin(thet)];
        end
    end
    
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

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
    
    if(ups==0&&sig>0)
        if(ext)
            activnet_pullext(N,tt,z0,zet,L,mu,kap,del,nu,psi,sig,D,Df,ncnt,lf,r,tinc,fileID);
        else
            activnet_pull(N,tt,z0,zet,L,mu,kap,del,nu,psi,sig,D,Df,ncnt,lf,r,tinc,fileID);
        end
    else
        activnet_act(N,tt,z0,zet,L,mu,kap,del,nu,psi,sig,D,Df,ncnt,lf,r,tinc,fileID);
    end
    
end
