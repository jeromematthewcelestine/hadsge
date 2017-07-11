function data = load_aggs_data()

output = gdpc1.Data(:,2);
consumption = pcecc96.Data(:,2);

[~, output_hp] = hpfilter(output);
[~, consumption_hp] = hpfilter(consumption);

% plot(output_trend); hold on; plot(output); hold off;
% plot(output_hp)
plot(output_hp); hold on; plot(consumption_hp); hold off;

data = [output_hp consumption_hp]';