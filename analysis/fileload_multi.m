%zet,L,mu,kap,lc,del,ups,phi,psi,r,sig,D,Df,ls,lf,tinc,tfin
%medium viscosity, filament length, extensional stiffness, bending
%stiffness, crosslinker spacing, crosslinker friction, myosin force, myosin
%fraction, spatial variation of ups, recycling rate, stress applied, domain
%size, portion of domainedge deciated to applying forces, 


r=0;
origbp = pwd;
%#ok<*ST2NM>
bp = '/Users/wmcfadden/activ_begin';
cd(bp);
files = dir;
files = {files.name};
stoA = [];
stoG0 = [];
stoB = [];
stoC = [];
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
    %         r=str2num(paree{10});sig=str2num(paree{10});D=str2num(paree{11});Df=str2num(paree{12});ncnt=str2num(paree{13});lf=str2num(paree{14});
            %         tinc=str2num(paree{16});tfin=str2num(paree{17});
%             sig = sig*sin(psi); % stupid line cause I messed something up temporarily
            fclose(fid);
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
                    code{1}
                    length(t)
                    if(length(t)>10)
                        lst = size(zt,1);
                        trp = repmat((1:lst)'/lst,1,3);
                        cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
                        h = figure('Position', [50, 100, 1200, 600]);
                        hold on
            %             cc = copper(size(zt,1));
                        whitebg('black')
                        set(gcf,'Color',[0 0 0])
                        set(gcf,'InvertHardcopy','off')
                        for ind = 1:ceil(size(zt,1)/10):size(zt,1)
                            p = reshape(zt(ind,:),[],2);
                            p = [mod(p(:,1),2*D),mod(p(:,2),D)];

                            netplot(p,L,lf,ls,D,cc(ind,:));
                            drawnow
                        end


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

%                         ss = [];
%                         tt = [];
%                         for k=1:size(xt,1)
%                             if(sum(isnan(xt(k,:)))==0)
%                                 sv = bindata(dy(k,:),xt(k,:),brng);
%                                 subs = sp>Df*D&sp<D*(1-Df);
%                                 tt=[tt; t(k)];
%                                 ss=[ss; -nanmean(sv(subs)./sp(subs))];
%                             end
%                         end
%                         a_t = 1;
%                         b_t = 1;
%                         as = [];
%                         bs = [];
%                         cs = [];
%                         gofs= [];
%                         confs = [];
%                         for i=1:10
%                             [fito, gof] = fit(tt((end/2):end),ss((end/2):end),'a*x^b+c','StartPoint', [a_t, b_t,1],'Lower',[0 0 -Inf],'Display','final');
%                             a_t = fito.a*(1+0.5*randn);
%                             b_t = fito.b*(1+0.5*randn);
%                             as = [as;fito.a];
%                             bs = [bs;fito.b];
%                             cs = [cs;fito.c];
%                             gofs=[gofs;gof.rmse];
% %                             conf = confint(fito);
% %                             confs = [confs; diff(conf)./coeffvalues(fito)];
%                         end
                        if(1)
%                             spt = find(gofs==min(gofs));
%                             spt = spt(1);
%                             a = as(spt);
%                             b = bs(spt);
%                             c = cs(spt);
%                             
%                             kr0 = lc.^2./zet./del./(L-2*lc).^2*16*3/pi;
%                             kr2 = mean(diff(ss((end/2):end))./diff(tt((end/2):end)))/sig;
%                             kr3 = mean(diff(ss((3*end/4):end))./diff(tt((3*end/4):end)))/sig;
%                             
% %                             stoA = [stoA; a];
% % %                             stoG0 = [stoG0; G0];
% %                             stoB = [stoB; b];
% %                             stoC = [stoC; c];
% % %                             stokr = [stokr; kr];
%                             stokr2 = [stokr2; kr2];
%                             stokr3 = [stokr3; kr3];
%                             stokr0 = [stokr0; kr0];
%                             stoall = [stoall; zet L mu kap lc del ups phi psi r sig D Df ls lf];
% %                             stoconf = [stoconf; confs(spt,:)];
% %                             stogofs = [stogofs; gofs(spt)];
%                             stogams = [stogams; ss(end)/tt(end)];
%                             stoname = [stoname; code{1}];

%                             axes('Position',[.75 .7 .12 .2])
%                             title('strain plot')
%                             box on
%                             plot(tt,ss,'yo');
% %                             plot(tt((end/2):end),ss((end/2):end),'yo');
% %                             hold on
% %                             plot(tt(1:(end/2)),ss(1:(end/2)),'go');
% %                             plot(tt,a*tt.^b+c,'m','LineWidth',2);
%                             set(gca,'fontsize',14)
                            h_leg=annotation('textbox', [0.75 0.2 0.12 0.45],'BackgroundColor',[1 1 1],...
                                'String',{code{1},['zeta = ' num2str(zet)],['L/lc = ' num2str(L/lc)],['mu = ' num2str(mu)],['kap = ' num2str(kap)],...
                                ['zet*del = ' num2str(zet*del)]});
%                             ['zet*del = ' num2str(zet*del)],['n_m/n_p = ' num2str(kr0./kr2)],['exp = ' num2str(b)]});
                            set(h_leg,'FontSize',16);
                            set(h,'PaperPositionMode','auto')
                            drawnow
                            print('-dpng','-r0',[code{1} '_fig.png']);%saveas(h,['fig_' code{2} '.png'],'png');
                        end
                        close(h)
%                         clear mov
%                         figure('Position', [50, 100, 1200, 600]);
%                         lst = size(zt,1);
%                         trp = repmat((1:lst)'/lst,1,3);
%                         cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
%                         temp=hot(2*lst);
%                         cc2 = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*temp(1:lst,:);
%                         indi = 1;
%                         for ind = 1:100:size(zt,1)
%                             p = reshape(zt(ind,:),[],2);
%                             p = [mod(p(:,1),2*D),mod(p(:,2),D)];
%                             whitebg('black')
%                             set(gcf,'Color',[0 0 0])
%                             set(gcf,'InvertHardcopy','off')
% 
%                             netplot_str(p,L,lf,ls,D,cc,cc2,0.25);
%                             mov(indi) = getframe;
%                             clf
%                             indi = indi +1;
%                         end
%                         movie2avi(mov,[code{1} '_mov.avi']);
%                         close(gcf);
                    end
                end
            end
        end
    end
%     clearvars -except origbp bp files 
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
% save('fitvals','stokr','stokr2','stoall','stoconf','stogofs','stogams','stoname');
% 
zet = stoall(:,1);
L = stoall(:,2);
mu = stoall(:,3);
kap = stoall(:,4);
lc = stoall(:,5);
del = stoall(:,6);
ups = stoall(:,7);
phi = stoall(:,8);
psi = stoall(:,9);
r = stoall(:,10);
sig = stoall(:,11);
D = stoall(:,12);
