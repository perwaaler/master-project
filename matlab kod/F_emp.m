function out = F_emp(x, data)
% empirical distribution function of sample data evaluated at x, which can
% be a single point or a vector
dimx = size(x);
dimdata = size(data);

if dimx(1)>dimx(2)
 x = x';
end
if dimdata(1)>dimdata(2)
 data = data';
end

nn = length(x); 

ordered_sample = sort(data);
ordered_sample = ordered_sample(ones(1,nn),:);
interval_index = sum(x'-ordered_sample >= 0, 2);

out = (interval_index)/length(data);
end
