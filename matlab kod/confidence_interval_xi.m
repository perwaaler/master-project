sigma = 0.2255732;
xi = -0.4040614;
u = -.5;
n=204;
negated_pot = -danger_CR(find(-danger_CR>u))-u;
% [nlogL,acov] = gplike([xi,sigma], negated_pot) % computes the 
if xi<0
    ue = u - sigma/xi;
else
    ue = inf;
end

sigma_double_prime = @(y, h) -( log(gppdf(y,xi,sigma+h,0)) - 2*log(gppdf(y,xi,sigma,0)) + log(gppdf(y,xi,sigma-h,0)) )/h^2
xi_double_prime = @(y, h) -( log(gppdf(y,xi+h,sigma,0)) - 2*log(gppdf(y,xi,sigma,0)) + log(gppdf(y,xi-h,sigma,0)) )/h^2

sigma2 = @(y,xi,sigma) (sigma^2 - 2*sigma*y - xi*y.^2)./( sigma^2 * (sigma + xi*y).^2 )   %1/sigma^2 - (1 + xi)*( 2*sigma + xi*(x - u) ).*(x - u)./(sigma^2 + sigma*xi*(x - u)).^2;
integrandsigma = @(x,xi,sigma,u) sigma2(x,xi,sigma,u).*gppdf(x,xi,sigma,u);
V11 = integral(@(x) integrandsigma(x,xi,sigma,u),u ,ue , 'AbsTol',1e-10)*n



xi2 = @(x,xi,sigma,u) -2/xi^3*log( 1+(x-u)/sigma * xi ) + (2*sigma./(x-u) + 3*xi + xi^3)./(xi^2*(sigma./(x-u) + xi).^2); % second derivative w.r.t. xi
integrandxi = @(x,xi,sigma,u) xi2(x,xi,sigma,u).*gppdf(x,xi,sigma,u);
V22 = -integral(@(x) integrandxi(x,xi,sigma,u),u,ue)*n

xisigma = @(x,xi,sigma,u) -(x-u).^2 * (1 + xi) ./ (xi*(sigma + (x-u)*xi).^2); %( log(1+xi/sigma*(x-u))/xi^2 - ( (x-u)-sigma )./(sigma^2 + sigma*xi*(x-u))    ).^2 %-(x-u).^2 * (1 + xi) ./ (xi*(sigma + (x-u)*xi).^2); %
integrandxisigma = @(x,xi,sigma,u) xisigma(x,xi,sigma,u).*gppdf(x,xi,sigma,u);
V12 = -integral(@(x) integrandxisigma(x,xi,sigma,u),u,ue)*n

V = [[V11,V12];[V12,V22]]

%V^-1