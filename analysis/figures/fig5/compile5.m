
bns = 51;
ll = 4;
rl = 16;

allt = [];
allp = [];
allg = [];
allf = [];
alle = [];
allc = [];
allfe = [];
allfc = [];
alla = [];
allw = [];
alln = {};

bp = '/Users/wmcfadden/extend_rec_sweep_a/';
cd(bp);
files = dir;
files = {files.name};

for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end

bp = '/Users/wmcfadden/extend_rec_sweep_b/';
cd(bp);
files = dir;
files = {files.name};

for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end

bp = '/Users/wmcfadden/extend_rec_sweep_c/';
cd(bp);
files = dir;
files = {files.name};

for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end

if(size(allt)>0)
    save('allmeas','allt','allp','allg','alla','allf','alle','allc','alln')
end