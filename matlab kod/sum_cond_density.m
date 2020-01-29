function out = sum_cond_density(x,beta)
out = zeros(1,length(x))
print(out)
for i=1:length(x)
    out(i) = sum(exppdf(x(i),1./beta))
end
end