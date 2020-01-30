function out = expSum(x,beta)
out = zeros(1,length(x));
for i=1:length(x)
    out(i) = sum(exppdf(x(i),1./beta));
end
end