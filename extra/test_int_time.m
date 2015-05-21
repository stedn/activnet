



lcs =[];
for ind = 1:5:size(zt,1)
    ind
    p = reshape(zt(ind,:),[],2);
    p = [mod(p(:,1),2*D),mod(p(:,2),D)];
    
    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);

    indL = 1:length(p);
    indL=indL(mod(indL,ncnt)~=0);
    
    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);

    subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)=subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)+2*D;
    subpL(subpL(:,2)<D/4&subpR(:,2)>3*D/4,2)=subpL(subpL(:,2)<D/4&subpR(:,2)>3*D/4,2)+D;
    subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)=subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)+2*D;
    subpR(subpR(:,2)<D/4&subpL(:,2)>3*D/4,2)=subpR(subpR(:,2)<D/4&subpL(:,2)>3*D/4,2)+D;
    
    % extend ends slightly for smooth falloff
    subv = subpR-subpL;
    subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
    subpL = subpL - l0*lf/2*subv;
    subpR = subpR + l0*lf/2*subv;
    
    
    XY = [subpL subpR];
    
    subXY = XY(:,1)>2*D|XY(:,2)>D|XY(:,3)>2*D|XY(:,4)>D;
    
    extXY = XY(subXY, :);
    
%     tsub = extXY(:,1)>2*D;
%     extXY(tsub,:)=extXY(tsub,:)-repmat([2*D 0 2*D 0],sum(tsub),1);
%     tsub = extXY(:,3)>2*D;
%     extXY(tsub,:)=extXY(tsub,:)-repmat([2*D 0 2*D 0],sum(tsub),1);
%     tsub = extXY(:,2)>D;
%     extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
%     tsub = extXY(:,4)>D;
%     extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    
    XY = [XY; extXY];
    
    indL = [indL indL(subXY)];
    tic
    g = lineSegmentGrid(indL,XY,D,l0);
    lcs = [lcs length(g)];
end

