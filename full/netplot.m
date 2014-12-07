function netplot(p,L,lf,ls,D,clr)    
    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);

    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);
    subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)=subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)+2*D;
    subpL(subpL(:,2)<D/4&subpR(:,2)>3*D/4,2)=subpL(subpL(:,2)<D/4&subpR(:,2)>3*D/4,2)+D;
    subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)=subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)+2*D;
    subpR(subpR(:,2)<D/4&subpL(:,2)>3*D/4,2)=subpR(subpR(:,2)<D/4&subpL(:,2)>3*D/4,2)+D;
    
    subv = subpR-subpL;
    subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
    subpL = subpL - l0*lf/2*subv;
    subpR = subpR + l0*lf/2*subv;
    
    
    XY = [subpL subpR];
    subXY = XY(:,1)>2*D|XY(:,2)>D|XY(:,3)>2*D|XY(:,4)>D;
    
    extXY = XY(subXY, :);
    
    tsub = extXY(:,1)>2*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([2*D 0 2*D 0],sum(tsub),1);
    tsub = extXY(:,3)>2*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([2*D 0 2*D 0],sum(tsub),1);
    tsub = extXY(:,2)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    tsub = extXY(:,4)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    
    XY = [XY; extXY];
    line([XY(:,1)';XY(:,3)'],[XY(:,2)';XY(:,4)'],'Color',clr,'LineWidth',2);
    xlim([0 2*D]);
    ylim([0 D]);
end