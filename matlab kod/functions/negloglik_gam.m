function out = negloglik_gam(gpdpar, gampar, u)
% takes 
ss = 0;
for i=1:size(gampar,2)
    ss = ss - integral(@(x) gampdf(-x,gampar(1,i),gampar(2,i)) .* log(gppdf(x-(-u),gpdpar(2),gpdpar(1),0))  , -u, 0 );
out = ss;
end
