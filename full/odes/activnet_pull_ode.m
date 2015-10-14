function dz = activnet_pull_ode(t,z,zet,L,mu,kap,del,nu,psi,sig,D,Df,ncnt,lf,ext)
    
    %% compute intrafilament forces    
    l0 = L/(ncnt-1);
    p = reshape(z,[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];
    dp = zeros(size(p));
    for n=1:ncnt:length(p)
        va_orth=[0 0];
        va = [0 0];
        la = 0;
        for i=0:ncnt-2
            vb = mydiff(p(n+i,:),p(n+i+1,:),D);
            lb = sqrt(vb*vb');
            gam = (lb-l0)/l0;
            f = mu*vb/lb*gam;
            if(mu<0)
                f = -f*(1+99*gam^100/(0.05^100+gam^100));
            end
            dp(n+i,:) = dp(n+i,:) + f;
            dp(n+i+1,:) = dp(n+i+1,:) - f;
            vb_orth = [-vb(2) vb(1)];
            if(i>0)
                if(va_orth*vb'>0)
                    va_orth = -va_orth;
                end
                if(vb_orth*va'<0)
                    vb_orth = -vb_orth;
                end
                tor = kap/l0^2*acos(max(min(va*vb'/la/lb,1),0));
                dp(n+i-1,:)=dp(n+i-1,:)+tor*va_orth/la;
                dp(n+i,:)=dp(n+i,:)-tor*va_orth/la;
                dp(n+i+1,:)=dp(n+i+1,:)+tor*vb_orth/lb;
                dp(n+i,:)=dp(n+i,:)-tor*vb_orth/lb;
            end
            va = vb;
            va_orth = vb_orth;
            la = lb;
        end
    end
    
    
    
    %% add external force at centerline and constrain edges
    if(psi>0)
        val = sig*sin(psi*t);
    elseif(psi<0)
        val = sig*round(mod(0.5+-psi*t,1)).*(round(mod(0.55+-psi*t/2,1))-0.5)*2;
    else
        val = sig;
    end
    
    subp = p(:,1)>D-Df*D&p(:,1)<D+Df*D;
    ff = D-abs(p(subp,1)-D)/(1-Df);
    if(ext)
        dp(subp,1)=dp(subp,1) - D*val.*ff/sum(ff);
    else
        dp(subp,2)=dp(subp,2) - D*val.*ff/sum(ff);
    end
    
    subp = p(:,1)<Df*D;
    dp(subp,:)=dp(subp,:).*repmat(2*p(subp,1)/Df/D-1,1,2);

    subp = p(:,1)>2*D-Df*D;
    dp(subp,:)=dp(subp,:).*repmat(2*abs(p(subp,1)-2*D)/Df/D-1,1,2);

    subp = p(:,1)<Df*D/2|p(:,1)>2*D-Df*D/2;
    dp(subp,:)=0;

    %% and bring it all home
    
    dz = reshape(dp,[],1);
    
    
end