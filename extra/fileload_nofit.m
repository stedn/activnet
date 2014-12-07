r=0;
origbp = pwd;
%#ok<*ST2NM>
bp = '/Users/wmcfadden/active_lite/';
cd(bp);
files = dir;
files = {files.name};

pars = {};
for f = files
    if(strfind(f{1},'_scr.txt') )
        code = strsplit(f{1},'_');
        fid = fopen(f{1});
        C = textscan(fid, '%s','delimiter', '\n');
        pare = strsplit(C{1}{8}, '>');
        paree = strsplit(pare{1}, ' ');
        zet=str2num(paree{2});L=str2num(paree{3});mu=str2num(paree{4});kap=str2num(paree{5});lc=str2num(paree{6}); 
        del=str2num(paree{7});ups=str2num(paree{8});phi=str2num(paree{9});psi=str2num(paree{10});
%         r=str2num(paree{10});sig=str2num(paree{11});D=str2num(paree{12});Df=str2num(paree{13});ls=str2num(paree{14});lf=str2num(paree{15});
        r=str2num(paree{11});sig=str2num(paree{12});D=str2num(paree{13});Df=str2num(paree{14});ls=str2num(paree{15});lf=str2num(paree{16});
%         sig=str2num(paree{10});D=str2num(paree{11});Df=str2num(paree{12});ls=str2num(paree{13});lf=str2num(paree{14});
        fclose(fid);
        if(r==0.1||r==0.001)
            A = importdata([code{1} '_out.txt']);
            t = A(:,1);
            zt = A(:,2:end);
            if(length(t)>1)
                lst = size(zt,1);
                trp = repmat((1:lst)'/lst,1,3);
                cc = (1-trp.^2).*(winter(lst)*0.75+0.25*spring(lst))+(trp.^2).*copper(lst);
                h = figure('Position', [50, 100, 1200, 600]);
                hold on
                clear mov
                whitebg('black')
                set(gcf,'Color',[0 0 0])
                set(gcf,'InvertHardcopy','off')
    %             cc = copper(size(zt,1));
                indi=1;
                for ind = 1:100:size(zt,1)
                    p = reshape(zt(ind,:),[],2);
                    p = [mod(p(:,1),2*D),mod(p(:,2),D)];

                    netplot(p,L,lf,ls,D,cc(ind,:));
                    drawnow
                    mov(indi) = getframe;
                    indi=indi+1;
%                     clf
                end
%                 if(exist('mov'))
%                     movie2avi(mov,[code{1} '_mov.avi']);
%                 end
                h_leg=annotation('textbox', [0.75 0.2 0.12 0.45],'BackgroundColor',[1 1 1],...
                    'String',{code{2},['L = ' num2str(L)],['mu = ' num2str(mu)],['kap = ' num2str(kap)],['lc = ' num2str(lc)],...
                    ['del = ' num2str(del)], ['ups = ' num2str(ups)],['phi = ' num2str(phi)],['r = ' num2str(r)],...
                    ['sig = ' num2str(sig)],...['G = ' num2str(G)],['G0 = ' num2str(G0)],['xi = ' num2str(xi)],['xi0 = ' num2str(xi0)],['tau = ' num2str(tau)],...
                    ['end time = ' num2str(t(end))]});
                set(h_leg,'FontSize',16);
                set(h,'PaperPositionMode','auto')
                drawnow
                print('-dpng','-r0',[code{1} '_fig.png']);%saveas(h,['fig_' code{2} '.png'],'png');
                close(h)
            end
        end
    end
end

cd(origbp);
