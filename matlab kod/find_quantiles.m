function out = find_quantiles(x, data, weights)
% function that calculates the empirical quantiles - evaluated at points in x - corresponding to the 
% empirical cdf's of data, where data is a matrix of observations, each row 
% corresponding to a sample. Weights is an argument that determines the weight
% given to each cdf as they get added together to form the empirical cdf.
% of the exceedence data. note that 0<x<1. NOte that data excesses only,
% and have nan's for elements where data did not exceed threshold.

% find the points where the cdf jumps
jump_points = sort(data(:),1)';
jump_points = jump_points(find(isnan(jump_points)==0 ));

%compute emirical cdf of excesses
femp = weights'*F_emp(jump_points, data);

interval_index = sum(x' >= femp , 2);

out = jump_points(interval_index+1)
end


% clf
% F_inv = @(y,scale,shape) scale/shape*((1-y).^-shape - 1)
% plot(F_inv(x, param(1),param(2)),jump_points(interval_index+1),'.')
% hold on
% plot(jump_points(interval_index+1),jump_points(interval_index+1)); hold off
% 
% 
% 
% plot(jump_points, femp)

