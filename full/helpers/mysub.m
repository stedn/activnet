function d = mysub(x1,x2,S)
    d = x2-x1;
    d(d>S/2) = d(d>S/2) - S;
    d(d<-S/2) = d(d<-S/2) + S;
end