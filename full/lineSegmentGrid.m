function g = lineSegmentGrid(indL,XY,D,l0)
persistent ga;

if isempty(ga)
   ga = zeros(10000000,5);
end

dg=D/floor(D/2/l0);
XYg = ceil(XY/dg);
% n_max = 100;
% grid = zeros(2*D/dg,D/dg,n_max);

xmin = min(XYg(:,1),XYg(:,3));
xmax = max(XYg(:,1),XYg(:,3));
ymin = min(XYg(:,2),XYg(:,4));
ymax = max(XYg(:,2),XYg(:,4));

str = 1;
for xn=1:2*D/dg
    for yn=1:D/dg
        subXY = xmin<=xn&xmax>=xn&ymin<=yn&ymax>=yn;
%         subXY = grid(xn,yn);
        gsub = lineSegmentIntersect(indL(subXY),XY(subXY,:));
        if(~isempty(gsub))
            lst = str+size(gsub,1)-1;
            ga(str:lst,:) = gsub;
            str = lst+1;
        end
    end
end

g = unique(ga(1:str-1,:),'rows');
end

