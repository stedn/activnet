bp = '/Users/wmcfadden/activ_free/';
cd(bp);
load('allmeas')




plot(allt(end,2:end)/10,(sqrt(allg(end,2:end))-sqrt(allg(end,2)))/sqrt(allg(end,2)))
hold on
plot(allt(end,2:end)/10,allc(end,2:end))
plot(allt(end,2:end)/10,alle(end,2:end))
