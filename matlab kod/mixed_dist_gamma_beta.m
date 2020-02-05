N_samples = 300;  % number of samples drawn from each sample
N = 100;         % number of "encounters"
gam_mean = 5;
gam_var = 0.01;
% parametes as fcn of mean and variance
a = gam_mean^2/gam_var;
b = gam_var/gam_mean;
% generate sample
alfa = 0.00001; % shift needed to avoid 0
gam_par = gamrnd(a*ones(1,N), b*ones(1,N)) + alfa;
beta = 5.5;
X = betarnd(gam_par' * ones(1, N_samples), beta);
% compute true value of estiamand
mu = 0.05; % we are estimating the probability that x < mu
integrand = @(s) betacdf(mu, s, beta).*gampdf(s - alfa, a, b)

p_mu = integral(integrand, alfa,Inf)
% hist(X(:),50)
%% transform data, and plot to find appropriate max and min thresholds
transinv = @(x) 1./(x+0.051).^4;
transneg = @(x) -x;
transex = @(x) exp(-10*x);
trans = @(x) transneg(x);
trans_data = trans(X);
u_l = 0.08;
u_u = 0.12;
u_l_trans = trans(u_u);
u_u_trans = trans(u_l);
clf; plot(max(trans_data'),'.'); hold on; plot(ones(1,N)*u_l_trans,'r') 
plot(ones(1,N)*u_u_trans,'r'); plot(ones(1,N)*trans(mu),'g')

%% find initial values
trans_data = trans_data(:);
exceed = trans_data(find(trans_data > u_l_trans));
negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
par_init = fminsearch(@(par) negL(par, exceed, u_l_trans), [3 0.5])
%%
trans_data = trans_data(:);
negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
logplot = 0;
m = 10;
U = sort(trans(linspace(u_l,u_u,m)));
param_save = zeros(2,m);
ue_save = zeros(1,m)*nan;
p_nea = zeros(1,m); 
max_data = max(max(trans_data));
for k=1:m
    exceed = trans_data(find(trans_data > U(k)));
    par_init_temp = par_init;
    pu = length(exceed)/length(trans_data)
    param = fminsearch(@(par) negL(par, exceed, u_l_trans), par_init_temp);
    %while param = par_init_temp
    param_save(:,k) = param';
    ue = U(k) - param(1)/param(2);
    if param(2) < 0
        ue_save(k) = ue;
    end
    p_nea(k) = pu*(1 - gpcdf(trans(mu), param(2), param(1), U(k)));
end

clf
subplot(221)
plot(param_save(2,:))
title('estimate of xi')
subplot(222)
plot(log(p_nea),'r'); hold on; plot(ones(1,m)*log(p_mu),'g'); 
plot(ones(1,m)*log(sum(trans_data(:)>trans(mu))/length(trans_data)),'b')
title('logarithm of estimates of pu')
subplot(223)
plot(ue_save); hold on; plot(ones(1,m)*max_data)
title('upper endpoint estimates')
subplot(224)
plot(param_save(1,:))
title('estimates of sigma')

