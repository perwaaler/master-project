function out = integrand(x,par,beta,u)
out = expSum(x, beta,u).*pareto(x,par,u);
end