bp = '/Users/wmcfadden/activ_rec/';
cd(bp);
files = dir;
files = {files.name};
allt = [];
allp = [];
allg = [];
alla = [];
allf = [];
alln = {};
for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measure4
            alln = {alln{:} code};
        end
    end
end
save('allmeas','allt','allp','allg','alla','allf','alln')
