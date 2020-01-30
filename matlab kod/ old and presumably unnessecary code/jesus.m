function out = sum_cond_density(x,beta)
    out = sum(exppdf(x,1./beta));
end