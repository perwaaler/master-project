% The purpose of this code is to show an example of a situation in which
% applying some kind of inversion to the data improves the quality of the
% estimate obtained when using POT method. The data is sampled from a beta
% distribution picked specifically so that the situation mimics that of the
% traffic-security situation; i.e. we are looking to estimate the
% probability of the data falling below a certain threshold, and the
% probability of that happening is extremely low, resulting in the probability
% of dipping below such a threshold often being estimated to be zero.


N = 500;           % number of iterations
n = 500;           % number of samples generated each iteration

% parameter values of the beta distribution
a=6;
b=5.5;

% initial values when minimizing nagtive log-likelihood fcn.
beta_init = [1 -1];
inv_init = [1 1];                                                                                                 % plot(linspace(0,1,1000), betapdf(linspace(0,1,1000),a,b) )

p = betacdf(0.08,a,b);  % true value of probability we are trying to estimate = 6.0178e-05
psave = zeros(2,N);     % saves estimates for each iteration
trans = @(x) 1./x.^4    % specifies the type of transformation to use

u0 = -0.08;             % we are estimating P(X>u0)
u0_inv = trans(-u0);    % corresponding threshold for transformed data
u_beta = -0.4;          % threshold used when applying POT method
u_inv = trans(-u_beta); % corresponding threshold for transformed data


for k=1:N
    beta_sample = betarnd(a*ones(1,n), b*ones(1,n));
    inv_sample = trans(beta_sample);  % inverted data sample
    beta_sample = -beta_sample;       % negate data so that we can apply POT method

    % find exceedences
    beta_exceed = beta_sample(find(beta_sample > u_beta));
    inv_exceed = inv_sample(find(inv_sample > u_inv));

    negL = @(par, data, u) -sum( log(gppdf(data,par(2),par(1),u)) );

    par_beta = fminsearch(@(par) negL(par, beta_exceed, u_beta), beta_init)
    par_inv = fminsearch(@(par) negL(par, inv_exceed, u_inv), inv_init)

    pu = length(beta_exceed)/n;
    p0_beta = pu*(max(0,1 + par_beta(2)*(u0 - u_beta)/par_beta(1)) )^(-1/par_beta(2))

    pu = length(inv_exceed)/n;
    p0_inv = pu*(max(0,1 + par_inv(2)*(u0_inv - u_inv)/par_inv(1)) )^(-1/par_inv(2))

    % save estimates
    psave(:,k) = [p0_beta;p0_inv];
end

%%
m_beta = mean(psave(1,:));
m_inv = mean(psave(2,:));
avg_error_beta = sqrt(sum((psave(1,:)-p).^2)/N)
avg_error_inv = sqrt(sum((psave(2,:)-p).^2)/N)
clf
subplot(211)
plot(psave(1,:),'.'); hold on; plot(ones(1,N)*p); plot(ones(1,N)*m_beta)
title('estimate using non-transformed data')
subplot(212)
plot(psave(2,:),'.')
title('estimate using transformed data'); hold on; plot(ones(1,N)*p); plot(ones(1,N)*m_inv)
