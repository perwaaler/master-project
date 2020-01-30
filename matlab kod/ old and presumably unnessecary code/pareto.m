function out = logpareto(x,par,u,xo)
% logarithm of paretdo pdf. x0 represents the threshold such that if xi is lower than
% this threshold the exponential distribution is used instead
sigma = abs(par(1));
xi = par(2);
if xi>xo
    ue = u - sigma/xi;
    if x<ue
        out = -log(sigma)-(1+xi)/xi*log(1 + xi/sigma*(x - u));
    else
        out = 0.1
end
if abs(xi) <= xo
    out = -log(sigma)-(x-u)/sigma;
end
end