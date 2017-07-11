function download_aggs_data()

clear;

c = fred();

gdpc1 = fetch(c, 'GDPC1');
pcecc96 = fetch(c, 'PCECC96');

clear c;
save('aggs_data.mat');

