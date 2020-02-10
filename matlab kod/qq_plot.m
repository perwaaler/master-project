function qq_plot(exceed_data,scale,shape,u, iteration)
F_inv = @(y,scale,shape) scale/shape*((1-y).^-shape - 1);
N =  length(exceed_data);
excesses = exceed_data - u;
plot(sort(excesses), F_inv( [1:N] /(N+1), scale,shape), '.')
hold on
plot(sort(excesses), sort(excesses)); hold off
title('qq-plot')
end
