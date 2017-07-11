function data = load_aggs_data()

if (~exist('aggs_data.mat','file'))
    download_aggs_data();
end

load('aggs_data.mat');

output = gdpc1.Data(:,2);
consumption = pcecc96.Data(:,2);

[~, output_hp] = hpfilter(output, 1600);
[~, consumption_hp] = hpfilter(consumption, 1600);

data = [consumption_hp output_hp]';