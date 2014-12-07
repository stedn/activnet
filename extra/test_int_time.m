


ncnt = 10;
lf = 0.0;

D = 10;

L = 1;
lcs =[];
tms =[];
Ns = floor([logspace(1.5,3.5,15) logspace(1.5,3.5,15) logspace(1.5,3.5,15) logspace(1.5,3.5,15) ]);
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
    
    tsub = extXY(:,1)>2*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([2*D 0 2*D 0],sum(tsub),1);
    tsub = extXY(:,3)>2*D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([2*D 0 2*D 0],sum(tsub),1);
    tsub = extXY(:,2)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    tsub = extXY(:,4)>D;
    extXY(tsub,:)=extXY(tsub,:)-repmat([0 D 0 D],sum(tsub),1);
    
    XY = [XY; extXY];
    
    indL = [indL indL(subXY)];
    tic
    g = lineSegmentWrap(indL,XY,D,l0);
    tms = [tms toc];
    lcs = [lcs length(g)];
end
figure(5)
plot(Ns,tms,'b.');
figure(6)
plot(Ns,lcs,'b.');
