function Mz = activnet_mass(t,z,zet,L,mu,kap,xi,nu,psi,sig,Dx,Dy,Df,Dw,ncnt,lf,ext)

    %% create velocity coupling matrix    
    l0 = L/(ncnt-1);

    p = reshape(z,[],2);
    p = [mod(p(:,1),Dx),mod(p(:,2),Dy)];
    
    Mo = zeros(length(p));
    
    
    indL = 1:length(p);
    indL=indL(mod(indL,ncnt)~=0);
    
    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);

    subpL(subpL(:,1)<Dx/4&subpR(:,1)>3*Dx/4,1)=subpL(subpL(:,1)<Dx/4&subpR(:,1)>3*Dx/4,1)+Dx;
    subpL(subpL(:,2)<Dy/4&subpR(:,2)>3*Dy/4,2)=subpL(subpL(:,2)<Dy/4&subpR(:,2)>3*Dy/4,2)+Dy;
    subpR(subpR(:,1)<Dx/4&subpL(:,1)>3*Dx/4,1)=subpR(subpR(:,1)<Dx/4&subpL(:,1)>3*Dx/4,1)+Dx;
    subpR(subpR(:,2)<Dy/4&subpL(:,2)>3*Dy/4,2)=subpR(subpR(:,2)<Dy/4&subpL(:,2)>3*Dy/4,2)+Dy;
    
    % extend ends slightly for smooth falloff
    subv = subpR-subpL;
    subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
    subpL = subpL - l0*lf/2*subv;
    subpR = subpR + l0*lf/2*subv;
    
    XY = [subpL subpR];
    
    subXY = XY(:,1)>Dx|XY(:,2)>Dy|XY(:,3)>Dx|XY(:,4)>Dy;
    
    extXY = XY(subXY, :);
    
    tsub = extXY(:,1)>Dx;
    extXY(tsub,:)=extXY(tsub,:)-repmat([Dx 0 Dx 0],sum(tsub),1);
    tsub = extXY(:,3)>Dx;
    extXY(tsub,:)=extXY(tsub,:)-repmat([Dx 0 Dx 0],sum(tsub),1);
    tsub = extXY(:,2)>Dy;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 Dy 0 Dy],sum(tsub),1);
    tsub = extXY(:,4)>Dy;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 Dy 0 Dy],sum(tsub),1);
    
    XY = [XY; extXY];
    
    indL = [indL indL(subXY)];
    
    g = lineSegmentGrid(indL,XY,Dx,Dy,l0);

    f = min(1,max(0,(g-lf/2)/(1-lf)));
    for ind=1:size(g,1)
        i = g(ind,3);
        j = g(ind,4);
            
        edg = 1;
        
        if(g(ind,1)<lf)
            edg = edg*g(ind,1)/lf;
        elseif((1-g(ind,1))<lf)
            edg = edg*(1-g(ind,1))/lf;
        end

        if(g(ind,2)<lf)
            edg = edg*g(ind,2)/lf;
        elseif((1-g(ind,2))<lf)
            edg = edg*(1-g(ind,2))/lf;
        end
        
        Mo(i:i+1,i:i+1) = Mo(i:i+1,i:i+1) + edg*xi*[[1-f(ind,1) f(ind,1)]*(1-f(ind,1)); [1-f(ind,1) f(ind,1)]*f(ind,1)]; 
        Mo(i:i+1,j:j+1) = Mo(i:i+1,j:j+1) - edg*xi*[[1-f(ind,2) f(ind,2)]*(1-f(ind,1)); [1-f(ind,2) f(ind,2)]*f(ind,1)]; 
        
    end
    
    
    %% generate medium viscosity with assymetry for filament orientation
    
    vt = mydiff(p(2:end,:),p(1:end-1,:),Dx,Dy);
    v = [0 0; (vt(1:end-1,:)+vt(2:end,:))/2; 0 0];
    v(1:ncnt:end,:)=vt(1:ncnt:end,:);
    v(ncnt:ncnt:end,:)=vt(ncnt-1:ncnt:end,:);
    den = sqrt(v(:,1).^2+v(:,2).^2);
    gx = L*zet/ncnt*(1+abs(v(:,2))./den);
    gy = L*zet/ncnt*(1+abs(v(:,1))./den);
    
    %% and bring it all home
    
    Mz = (blkdiag(Mo+diag(gx),Mo+diag(gy)));
    
end