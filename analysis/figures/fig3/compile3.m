bp = '/Users/wmcfadden/extend_llc_ver/';
cd(bp);
files = dir;
files = {files.name};
allt = [];
allp = [];
allg = [];
allf = [];
alle = [];
allc = [];
alla = [];
allw = [];
alln = {};
for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            measureall
        end
    end
end
save('allmeas','allt','allp','allg','alla','allf','alle','allc','allw','alln')