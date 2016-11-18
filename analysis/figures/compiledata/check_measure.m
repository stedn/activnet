cd(bp)

%% load param file and decipher params
fid = fopen([bp code '_scr.txt']);
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
pare = strsplit(C{1}{9}, '>');
paree = strsplit(pare{1}, ' ');
paree = {paree{2:end}};
zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
xi=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
r=str2num(paree{11});sig=str2num(paree{12});Dx=str2num(paree{13});Dy=str2num(paree{14});Df=str2num(paree{15});
Dw=str2num(paree{16});ls=str2num(paree{17});lf=str2num(paree{18});nonl=str2num(paree{21});


%% load simulation data
t = [];
A = [];
if(mu<0)
    A = importdata([bp code '_out.txt']);
    if(isstruct(A))
        A = A.data;
        if(size(A,1)==1)
            imp2 = importdata([bp code '_out.txt'],' ',9);
            if(isfield(imp2,'data'))
                A = [A;imp2.data];
            end
        end
        t = A(:,1);
        zt = A(:,2:end);
   end
end

ups
length(t)