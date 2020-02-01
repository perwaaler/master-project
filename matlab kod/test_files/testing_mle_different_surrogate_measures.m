negL = @(par, data, u) -sum(log(gppdf(data,par(2),par(1),u)));
data = betarnd(ones(1,100)*1, ones(1,100)*1)
init = [1 0];
options = optimset('Display','iter','PlotFcns',@optimplotfval);
options.MaxIter = 100000;
options.MaxFunEval = 100000;
[par_beta, fcnval, exitflag] = fminsearch(@(par) negL(par, [data, 0), init, options)