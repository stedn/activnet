function fileload_cluster(bp)
    %zet,L,mu,kap,lc,del,ups,phi,psi,r,sig,D,Df,ls,lf,tinc,tfin
    %medium viscosity, filament length, extensional stiffness, bending
    %stiffness, crosslinker spacing, crosslinker friction, myosin force, myosin
    %fraction, spatial variation of ups, recycling rate, stress applied, domain
    %size, portion of domainedge deciated to applying forces, 
    origbp = pwd;
    %#ok<*ST2NM>
    cd(bp);
    files = dir;
    files = {files.name};
    stoG = [];
    stoG0 = [];
    stoxi = [];
    stoxi0 = [];
    stokr = [];
    stokr2 = [];
    stokr3 = [];
    stokre = [];
    stokr0 = [];
    stoall = [];
    stoconf = [];
    stogofs = [];
    stogams = [];
    stoname = [];


    pars = {};
    for f = files
        if(strfind(f{1},'_scr') )
            code = strsplit(f{1},'_');
            if(exist([code{1} '_out.txt'],'file'))
                fid = fopen(f{1});
                C = textscan(fid, '%s','delimiter', '\n');
                try
                    pare = strsplit(C{1}{9}, '>');
                catch
                    pare = strsplit(C{1}{8}, '>');
                end
                paree = strsplit(pare{1}, ' ');
                paree = {paree{2:end}};
                zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
                del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
                r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});
       
                
                sig = sig*sin(psi); % stupid line cause I messed something up temporarily
                fclose(fid);
                sig.*lc^2/zet/del/L^2
                if(1)
                    imp = importdata([code{1} '_out.txt'],' ',4);
                    A = imp.data;
                    if(~isempty(A))
                        if(size(A,1)==1)
                            imp2 = importdata([code{1} '_out.txt'],' ',9);
                            if(isfield(imp2,'data'))
                                A = [A;imp2.data];
                            end
                        end
                        t = A(:,1);
                        zt = A(:,2:end);
                        length(t)
                        if(length(t)>1)
                            lst = size(zt,1);


                            xt = zt(:,1:end/2);
                            yt = zt(:,end/2+1:end);
                %             dy = yt-repmat(yt(1,:),size(yt,1),1);

                            coff = abs(2*median(yt(2,:)-yt(1,:)));
                            dy = 0*yt;
                            for k=2:size(yt,1)
                                dyc = yt(k,:)-yt(k-1,:);
                                subs = abs(dyc)<coff;
                                dy(k,subs)=dy(k-1,subs)+dyc(subs);
                                dy(k,~subs)=dy(k-1,~subs);
                            end

                            nbins = 30;
                            brng = linspace(0,2*D,nbins+1)';
                            sp = brng(1:end-1)+brng(2)/2;
                            sp(sp>D)=2*D-sp(sp>D);

                            ss = [];
                            tt = [];
                            for k=1:size(xt,1)
                                if(sum(isnan(xt(k,:)))==0)
                                    sv = bindata(dy(k,:),xt(k,:),brng);
                                    subs = sp>Df*D&sp<D*(1-Df);
                                    tt=[tt; t(k)];
                                    ss=[ss; -nanmean(sv(subs)./sp(subs))];
                                end
                            end
                            G0 = 3*pi/8*mu/L*(L/lc + 2*lc/L - 3);
                            xi0 = zet*D^2/lc;
                            kr0 = lc^2/zet/del/L^2;
                            a_t = sig/G0;
                            b_t = G0/xi0;
                            c_t = G0*kr0;
                            as = [];
                            bs = [];
                            cs = [];
                            gofs= [];
                            confs = [];
%                             for i=1:10
%                                 [fito, gof] = fit(tt,ss,'a*(1-exp(-b*x)+c*x)','StartPoint', [a_t, b_t, c_t],'Lower',[0 0 0],'Display','final');
%                                 a_t = fito.a*(1+0.5*randn);
%                                 b_t = fito.b*(1+0.5*randn);
%                                 c_t = fito.c*(1+0.5*randn);
%                                 as = [as;fito.a];
%                                 bs = [bs;fito.b];
%                                 cs = [cs;fito.c];
%                                 gofs=[gofs;gof.rmse];
%                                 conf = confint(fito);
%                                 confs = [confs; diff(conf)./coeffvalues(fito)];
%                             end
                            if(~isnan(gofs))
                                spt = find(gofs==min(gofs));
                                spt = spt(1);
                                G = sig/as(spt);
                                xi = G/bs(spt);
                                kr = cs(spt)/G;
                                kr2 = mean(diff(ss((end/2):end))./diff(tt((end/2):end)))/sig;
                                kr3 = mean(diff(ss((3*end/4):end))./diff(tt((3*end/4):end)))/sig;

                                stoG = [stoG; G];
                                stoG0 = [stoG0; G0];
                                stoxi = [stoxi; xi];
                                stoxi0 = [stoxi0; xi0];
                                stokr = [stokr; kr];
                                stokr2 = [stokr2; kr2];
                                stokr3 = [stokr3; kr3];
                                stokr0 = [stokr0; kr0];
                                stoall = [stoall; zet L mu kap lc del ups phi psi r sig D Df ls lf];
                                stoconf = [stoconf; confs(spt,:)];
                                stogofs = [stogofs; gofs(spt)];
                                stogams = [stogams; ss(end)/tt(end)];
                                stoname = [stoname; code{1}];


                            end
                        end
                    end
                end
            end
        end
    end
    save('fitvals','stoG','stoxi','stokr','stokr2','stoall','stoconf','stogofs','stogams','stoname');
end
