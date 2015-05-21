function netplot_str(p,L,lf,ls,D,clr,clr2,str_max)    
    ncnt = ceil(L/ls)+1;
    l0 = L/(ncnt-1);

    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);
    subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)=subpL(subpL(:,1)<D/2&subpR(:,1)>3*D/2,1)+2*D;
    subpL(subpL(:,2)<D/3&subpR(:,2)>2*D/3,2)=subpL(subpL(:,2)<D/3&subpR(:,2)>2*D/3,2)+D;
    subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)=subpR(subpR(:,1)<D/2&subpL(:,1)>3*D/2,1)+2*D;
    subpR(subpR(:,2)<D/3&subpL(:,2)>2*D/3,2)=subpR(subpR(:,2)<D/3&subpL(:,2)>2*D/3,2)+D;
    
    
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
    l = l0 + l0*lf;
    str = sqrt((XY(:,3)-XY(:,1)).^2 + (XY(:,4)-XY(:,2)).^2);
    kn = (str-l)>=0;
    str = min(length(clr),ceil(length(clr)*abs(str-l)/l/str_max));
    for i =1:length(XY)
        if(kn(i))
            line([XY(i,1)';XY(i,3)'],[XY(i,2)';XY(i,4)'],'Color',clr(str(i),:),'LineWidth',1);
        else
            line([XY(i,1)';XY(i,3)'],[XY(i,2)';XY(i,4)'],'Color',clr2(str(i),:),'LineWidth',1);
        end
    end
    xlim([0 2*D]);
    ylim([0 D]);
end