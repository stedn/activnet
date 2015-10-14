function fileload_oscil_cluster(bp)
    %zet,L,mu,kap,lc,del,ups,phi,psi,r,sig,D,Df,ls,lf,tinc,tfin
    %medium viscosity, filament length, extensional stiffness, bending
    %stiffness, crosslinker spacing, crosslinker friction, myosin force, myosin
    %fraction, spatial variation of ups, recycling rate, stress applied, domain
    %size, portion of domainedge deciated to applying forces, 
    r=0;
    origbp = pwd;
    %#ok<*ST2NM>
    cd(bp);
    files = dir;
    files = {files.name};
    stoGp = [];
    stoGpp = [];
    stomega = [];
    stoall = [];
    stoconf = [];
    stogofs = [];
    stogams = [];



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
        %         r=str2num(paree{10});sig=str2num(paree{10});D=str2num(paree{11});Df=str2num(paree{12});ncnt=str2num(paree{13});lf=str2num(paree{14});
                %         tinc=str2num(paree{16});tfin=str2num(paree{17});
                fclose(fid);
                if(psi<1)
                    imp = importdata([code{1} '_out.txt'],' ');
                    A = imp.data;
                    if(~isempty(A))
                        t = A(:,1);
                        zt = A(:,2:end);
                        code{1}
%                         length(t)
                        if(length(t)>1)
                            lst = size(zt,1);
                            trp = repmat((1:lst)'/lst,1,3);
%                             cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
%                             h = figure('Position', [50, 100, 1200, 600]);
%     %                         hold on
%                 %             cc = copper(size(zt,1));
%                             whitebg('black')
%                             set(gcf,'Color',[0 0 0])
%                             set(gcf,'InvertHardcopy','off')
%                             for ind = 1:10000:size(zt,1)
%                                 p = reshape(zt(ind,:),[],2);
%                                 p = [mod(p(:,1),2*D),mod(p(:,2),D)];
% 
%                                 netplot(p,L,lf,ls,D,cc(ind,:));
%                                 drawnow
%                                 clf
%                             end


                            xt = zt(:,1:end/2);
                            yt = zt(:,end/2+1:end);
                %             dy = yt-repmat(yt(1,:),size(yt,1),1);

                            coff = abs(3*median(median(diff(yt))));
                            dy = 0*yt;
                            for k=2:size(yt,1)
                                dyc = yt(k,:)-yt(k-1,:);
                                subs = abs(dyc)<Inf;
                                dy(k,subs)=dy(k-1,subs)+dyc(subs);
                                dy(k,~subs)=dy(k-1,~subs);
                            end

                            nbins = 30;
                            brng = linspace(0,2*D,nbins+1)';
                            sp = brng(1:end-1)+brng(2)/2;
                            sp(sp>D)=2*D-sp(sp>D);

                            ss = [];
                            tt = [];
                            for k=1:10000%size(xt,1)
                                sv = bindata(dy(k,:),xt(k,:),brng);
                                subs = sp>Df*D&sp<D*(1-Df);
                                tt=[tt; t(k)];
                                ss=[ss; -nanmean(sv(subs)./sp(subs))];
                            end
                            G0 = 3*pi/8*mu/L*(L/lc + 2*lc/L - 3);
                            a_t = 1000*sig/G0;
                            b_t = psi;
                            c_t = 0;
                            as = [];
                            bs = [];
                            cs = [];
                            ds = [];
                            es = [];
                            gofs= [];
                            confs = [];
                            for i=1:10
                                [fito, gof] = fit(tt,ss,'a*sin(b*x+c)+d*x + e','StartPoint', [a_t, b_t, c_t, 0, 0],'Display','final');
                                a_t = fito.a*(1+0.5*randn);
                                b_t = fito.b*(1+0.5*randn);
                                c_t = fito.c*(1+0.5*randn);
                                as = [as;fito.a];
                                bs = [bs;fito.b];
                                cs = [cs;fito.c];
                                ds = [ds;fito.d];
                                es = [es;fito.e];
                                gofs=[gofs;gof.rmse];
                                conf = confint(fito);
                                confs = [confs; diff(conf)./coeffvalues(fito)];
                            end
    %                        
                            if(~isnan(gofs))
                                spt = find(gofs==min(gofs));
                                spt = spt(1);
                                Gp = sig/as(spt)*cos(cs(spt));
                                Gpp = sig/as(spt)*sin(cs(spt));

                                stoGp = [stoGp; Gp];
                                stoGpp = [stoGpp; Gpp];
                                stomega = [stomega; bs(spt)];
                                stoall = [stoall; zet L mu kap lc del ups phi psi r sig D Df ls lf];
                                stoconf = [stoconf; confs(spt,:)];
                                stogofs = [stogofs; gofs(spt)];
                                stogams = [stogams; ss(end)];
fprintf('%s\n','stoGp')
    fprintf('%f\n', stoGp)
    fprintf('%s\n','stoGpp')
    fprintf('%f\n', stoGpp)
    fprintf('%s\n','stomega')
    fprintf('%f\n', stomega)
    fprintf('%s\n','stoall')
    fprintf('%f\n', stoall)
    fprintf('%s\n','stoconf')
    fprintf('%f\n', stoconf)
    fprintf('%s\n','stogofs')
    fprintf('%f\n', stogofs)
    fprintf('%s\n','stogams')
    fprintf('%f\n', stogams)
    %                             axes('Position',[.75 .7 .12 .2])
    %                             title('strain plot')
    %                             box on
    %                             plot(tt,ss,'yo');
    %                             hold on
    %                             plot(tt,as(spt)*sin(bs(spt)*tt + cs(spt))+ds(spt)*tt + es(spt),'m','LineWidth',2);
    %                             plot(tt,sig/Gp*sin(psi*tt),'b','LineWidth',2);
    %                             set(gca,'fontsize',14)
    %                             h_leg=annotation('textbox', [0.75 0.2 0.12 0.45],'BackgroundColor',[1 1 1],...
    %                                 'String',{code{1},['zeta = ' num2str(zet)],['L = ' num2str(L)],['lc = ' num2str(lc)],['mu = ' num2str(mu)],['kap = ' num2str(kap)],...
    %                                 ['del = ' num2str(del)]});
    %                             set(h_leg,'FontSize',16);
    %                             set(h,'PaperPositionMode','auto')
    %                             drawnow
    %                             print('-dpng','-r0',[code{1} '_fig.png']);%saveas(h,['fig_' code{2} '.png'],'png');
                            end
%                             close(h)
                        end
                    end
                end
            end
        end
    end
    % minds = stogofs>0.1;
    % stoG(minds,:)=[];
    % stoG0(minds,:)=[];
    % stoxi(minds,:)=[];
    % stoxi0(minds,:)=[];
    % stokr(minds,:)=[];
    % stokr0(minds,:)=[];
    % stoall(minds,:)=[];
    % stoconf(minds,:)=[];
    % stogofs(minds,:)=[];
%     save('fitvals','stoGp','stoGpp','stomega','stoall','stoconf','stogofs','stogams');
    
end
% zet = stoall(:,1);
% L = stoall(:,2);
% mu = stoall(:,3);
% kap = stoall(:,4);
% lc = stoall(:,5);
% del = stoall(:,6);
% ups = stoall(:,7);
% phi = stoall(:,8);
% psi = stoall(:,9);
% r = stoall(:,10);
% sig = stoall(:,11);
% D = stoall(:,12);
% figure
% hold on
% for dm = [10000 100000 1000000 10000000]
%     
%         subind = del./mu==dm&L./lc==12.5;
%         X = log(psi(subind));
%         Y = log(abs(stoGp(subind))./mu(subind));
%         YY = log(abs(stoGpp(subind))./mu(subind));
%         plot(X,Y,'o');
%         plot(X,YY,'x');
%         p = polyfit(X,Y,2);
%         x1 = linspace(min(X),max(X),100);
%         y1 = polyval(p,x1);
%         plot(x1,y1)
%         p = polyfit(X,YY,2);
%         x1 = linspace(min(X),max(X),100);
%         y1 = polyval(p,x1);
%         plot(x1,y1,'--')
% end