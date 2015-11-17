bp = '/Users/wmcfadden/activ_please/';
cd(bp);
files = dir;
files = {files.name};
alltt = [];
alltest = [];
allstrain = [];
allstrun = [];
allvisc = [];
alldisp=[];
allenrg=[];
allangs=[];
details = [];
for f = files
    if(strfind(f{1},'_scr') )
        code = strsplit(f{1},'_');
        if(exist([code{1} '_out.txt'],'file'))
            code = code{1}
            fileload_one
        end
    end
end