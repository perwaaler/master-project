function out = p_exceed(u,s,data)
% computes the probability to exceed threshold u.
N = max(data.ID);
trans_data = exp(-s*data.TTCFEA);
sum_negloglik = 0;
P = zeros(1,N);
for k = 1:N
    index = find(data.ID == k);
    nn = length(index);                                                  % number of possible outcomes at moment of first evasive action
    trans_data_k = trans_data(index);
    trans_data_k = trans_data_k( find(trans_data_k>u) );
    if nn > 0
        P(k) = length(trans_data_k)/nn;
    end
end
out = sum(P)/N;
end