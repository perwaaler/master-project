lamda = [1 2 3]
exppdf(1,1./lamda)
con_density = @(x,beta) exppdf(x,1./beta)
%sum_cond_density = @(x,beta) sum(exppdf(x,1./beta))
%tit = @(x,par) parretopdf(x,par)*2

con_density(1,[1,2,3,4,5,6])
sum_cond_density(1,[1,2,1,4,5,6])



fun = @(x,c) 1./(x.^3-2*x-c);
q = integral(@(x)fun(x,5),0,2)


%% example of using fminsearch
fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1.2,1];
options = optimset('PlotFcns',@optimplotfval);
x = fminsearch(fun,x0,options)
%%
beta = betarnd(3,1,[1,1000]);
u = 10;
xo = 10^-8;
options = optimset('PlotFcns',@optimplotfval);
xinit = [2       1e-07]
[xx,fval]=fminsearch(@(par)NLL(par,beta,u,xo),xinit)









