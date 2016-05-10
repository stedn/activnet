bp = '/Users/wmcfadden/activ_domain_ed_2/';
cd(bp);
files = dir;
files = {files.name};
allt = [];
allp = [];
allg = [];
alla = [];
allc = [];
alle = [];
allf = [];
alln = {};
for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measure5
        end
    end
end
if(size(allt)>0)
    save('allmeas','allt','allp','allg','alla','allf','alle','allc','alln')
end