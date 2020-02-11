% code for computing and plotting empirical vs model 1uantiles

k=7
pp = (1:100)/101;                                                   % select percentiles
data_string = trans_data(:);                                          % turn the data into a single data vector
ordered_excesses = sort(data_string(find(data_string>u))) - u;        %   
u = U(k);                                                             % set desired threshold

% gather arguments for find_quantiles(x, data, weights)
exceed_data = trans_data(    find(max(trans_data,[],2) > u) ,:   );   % matrix with data ONLY from encounters with atleast one exceedence of threshold
pu_i = sum(exceed_data > u, 2)/ size(trans_data,2);                   % probability to exceed threshold for each encounter
weights = pu_i/sum(pu_i) ;                                            % weights for each distribution function.

% replace obs. below u with nan's
below_u = find(exceed_data <= u);
exceed_data(below_u) = nan*ones(1,length(below_u));

exceed_data = exceed_data - u;

% estimate model parameters
trans_data_copy = data_string;
exceed = trans_data_copy(find(trans_data_copy(:) > u_l_trans));
negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
param = fminsearch(@(par) negL(par, exceed, u), [param_save(1,k),param_save(2,k)]);

% inverse cdf of estimated gppd
F_inv = @(y,scale,shape) scale/shape*((1-y).^-shape - 1);

% plot empirical quaniles vs model quantiles
clf
plot(F_inv(pp,param(1),param(2)), find_quantiles(pp, exceed_data, weights),'.')
hold on
plot(F_inv(pp,param(1),param(2)), F_inv(pp,param(1),param(2))); hold off
title('qq-plot')