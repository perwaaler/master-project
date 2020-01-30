function out = negloglik_exp(param, u, s, data)
% computes negative loglikelihood function for exponentially transformed data a given threshold u.
% threshold). Param is a vector containing scale and shape parameter of the
% generalized pareto distribution. s is the power of the exponential
% transformation; exp(-s*x). data is a  table containing id of encounters
% in one column and potential ttc values in another.

N = max(data.ID);
trans_data = exp(-s*data.TTCFEA);
sum_negloglik = 0;
for k = 1:N
    index = find(data.ID == k);
    nn = length(index);                                                  % number of possible outcomes at moment of first evasive action
    trans_data_k = trans_data(index);
    trans_data_k = trans_data_k( find(trans_data_k>u) );
    term_k = -log(gppdf(trans_data_k, param(2), param(1), u));
    term_k = term_k( find(~isnan(term_k)) );
    if length(term_k) > 0
        sum_negloglik = sum_negloglik + sum(term_k)/nn;
    end
end
out = sum_negloglik;
end