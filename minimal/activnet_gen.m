%% setup initiation params
function activnet_gen(zet,L,mu,kap,lc,del,ups,phi,psi,r,sig,D,Df,ls,lf,tinc,tfin)

    %% convert inputs
    ncnt = ceil(L/ls)+1;
    N = floor(4*D^2/lc/L);
    
    nu=[];
    if(ups>0)
        nu = ups*double(rand(N,N)<phi);
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
    
    activnet_act(N,tt,z0,zet,L,mu,kap,del,nu,psi,sig,D,Df,ncnt,lf,r,tinc,fileID);
    
end