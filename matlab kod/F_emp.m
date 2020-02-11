function out = F_emp(x, data)
% empirical distribution function of sample data evaluated at x. x can be a
% point or a row vector, and data can be a row vector or a matrix where
% each row contains a sample from some distribution. Returns a matrix where
% row i is the empirical d.f. corresponding to row_i(data) evaluated at the
% points in x. if the samples contains different number of obeservations,
% sample matrix should be padded with nan's.

sample_size = size(data,2);          % size of each sample, nan's included
nan_ind = isnan(data);               % index of nan's
n_nan = sum(nan_ind,2);              % number of nan's for each row
data(nan_ind) = zeros(1,sum(n_nan)); % replace nan's with zeros

if size(x,1)>size(x,2)
 x = x';
end

n_eval = size(x,2);          % number of points for which the emp. dist. gets evaluated
n_samples = size(data,1);    % number of samples

% correction term and factor
nan_sub = (n_nan./sample_size * ones(1,n_eval))';
nan_sub = nan_sub(:);
nan_mult = 1./(1 - nan_sub); 

n_copies = ((1:n_samples)'*ones(1,n_eval))'; 
n_copies = n_copies(:);

xx = x(ones(1,n_samples),:)';
xx = xx(:);
ordered_sample = sort(data,2);
ordered_sample = ordered_sample(n_copies,:); % constructs a matrix where the first nn rows are row 1, followed by nn copies of row 2 etc...
interval_index = sum(xx >= ordered_sample , 2);

out = (interval_index/sample_size - nan_sub).*nan_mult;
out = vec2mat(out,n_eval);


end
