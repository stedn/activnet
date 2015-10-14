


ncnt = 5;
lf = 0.0;

D = 4;

L = 1;
lcs =[];
Ns = floor([logspace(1.5,2.5,15) logspace(1.5,2.5,15) logspace(1.5,2.5,15) logspace(1.5,2.5,15) ]);
for N=Ns
N
    p = zeros(N*ncnt,2);
    for i=1:N
        p((i-1)*ncnt+1,:) = D*[rand rand];
        thet = rand*2*pi;
        for j = 2:ncnt
            p((i-1)*ncnt+j,:) = p((i-1)*ncnt+j-1,:)+L/(ncnt-1)*[cos(thet) sin(thet)];
            thet = thet+pi/16*randn;
        end

    end
    
    l0 = L/(ncnt-1);

    subpL = p;
    subpR = p;
    subpL=subpL(mod(1:length(subpL),ncnt)~=0,:);
    subpR=subpR(mod(1:length(subpR),ncnt)~=1,:);
    subpL(subpL(:,1)<D/4&subpR(:,1)>3*D/4,1)=subpL(subpL(:,1)<D/4&subpR(:,1)>3*D/4,1)+1*D;
    subpL(subpL(:,2)<D/4&subpR(:,2)>3*D/4,2)=subpL(subpL(:,2)<D/4&subpR(:,2)>3*D/4,2)+D;
    subpR(subpR(:,1)<D/4&subpL(:,1)>3*D/4,1)=subpR(subpR(:,1)<D/4&subpL(:,1)>3*D/4,1)+1*D;
    subpR(subpR(:,2)<D/4&subpL(:,2)>3*D/4,2)=subpR(subpR(:,2)<D/4&subpL(:,2)>3*D/4,2)+D;
    
    subv = subpR-subpL;
    subv = subv./repmat(sqrt(subv(:,1).^2+subv(:,2).^2),1,2);
    subpL = subpL - l0*lf/2*subv;
    subpR = subpR + l0*lf/2*subv;
    
    
    XY = [subpL subpR];
    subXY = XY(:,1)>1*D|XY(:,2)>D|XY(:,3)>1*D|XY(:,4)>D;
    
    extXY = XY(subXY, :);
    
    tsub = extXY(:,1)>1*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([1*D 0 1*D 0],sum(tsub),1);
    tsub = extXY(:,3)>1*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([1*D 0 1*D 0],sum(tsub),1);
    tsub = extXY(:,2)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    tsub = extXY(:,4)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    
    XY = [XY; extXY];
    
    subXY = XY(:,1)>1*D|XY(:,2)>D|XY(:,3)>1*D|XY(:,4)>D;

    extXY = XY(subXY, :);

    tsub = extXY(:,1)>1*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([1*D 0 1*D 0],sum(tsub),1);
    tsub = extXY(:,3)>1*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([1*D 0 1*D 0],sum(tsub),1);
    tsub = extXY(:,2)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    tsub = extXY(:,4)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);

    XY = [XY; extXY];
    
    out = lineSegmentIntersect(XY,XY);

    g = out.intNormalizedDistance1To2;
    f = min(1,max(0,(g-lf/2)/(1-lf)));
    [ns, ms] = find(out.intAdjacencyMatrix);

    ltot = length(ns);
    lcs = [lcs L*N/ltot];
end
plot(Ns,lcs,'g.');