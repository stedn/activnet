function v = mydiff(p1,p2,D)
    v = [mysub(p1(:,1),p2(:,1),2*D) mysub(p1(:,2),p2(:,2),D)];
end