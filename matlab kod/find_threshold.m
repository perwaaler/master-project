function out = find_threshold(data, lower_frac, upper_frac)
% function that computes the lower and upper threshold based on what
% fraction of the data points should be used in the analysis, so that
% approximately lower_frac and upper_frac of the data falls above lower_u
% and upper_u respectively.


sorted_data = sort(data);
lower_u = sorted_data(floor( lower_frac*length(data) ) );
upper_u = sorted_data(floor( upper_frac*length(data) ) );

out = [lower_u, upper_u];
end