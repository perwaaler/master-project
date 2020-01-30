function out = NLL(par, beta, u,xo)
sigma = abs(par(1));
xi = par(2);
integrand = @(x,par,beta,u,xo) expSum(x, beta).*pareto(x,par,u,xo);
if xi > -xo
    out = -integral(@(x)integrand(x,par,beta,u,xo), u, Inf);
end
if xi <= -xo
    ue = u - sigma/xi;
    out = -integral(@(x)integrand(x,par,beta,u,xo), u, ue);
end
end