
% The purpose of this code is to experiment with a known distribution to
% see if applying a transformation to the changes the estimate of the
% probability to exceed a certain threshold using EVT. More specifically,
% the distribution is U(0,1), and the transformation is picked specifically
% such that the resulting data is exp(1) distributed, and thus falls into
% the GEV family, making the assumption that the exceedences are GPPD
% distributed precise.

N = 10000;           % number of outer loops
n = 500;            % number of datapoints generated each loop
psave = zeros(2,N); % saves estimates of probabily for each method


u_uni = 0.8;                % threshold used when performing EVT analysis
u_exp = -log(1 - u_uni);    % corresponding threshold for transformed data
u0 = 0.98;                  % we are estimating P(X>u0). Note that the true probability is 1 - u0.
u0_trans = log(1/(1 - u0)); % corresponding value for transformed data

% initial values for log-likelihood minimization
exp_init = [1 0]; 
uni_init = [1 -0.1];

for k = 1:N
    uni_sample = rand(1,n);
    exp_sample = -log(1 - uni_sample);
    
    % find exceedences
    uni_exceed = uni_sample(find(uni_sample>u_uni));
    exp_exceed = exp_sample(find(exp_sample>u_exp));


    negL = @(par, data, u) -sum( log(gppdf(data,par(2),par(1),u)) );
    uni_init_temp = uni_init;
    exp_init_temp = exp_init;
    par_uni = fminsearch(@(par) negL(par, uni_exceed, u_uni), uni_init_temp);
    par_exp = fminsearch(@(par) negL(par, exp_exceed, u_exp), exp_init_temp);
    
    % try new initial guess in case of bad initial guess
    while par_uni == uni_init_temp
        uni_init_temp = [max(0.1, uni_init(1) + normrnd(0, 1.4^2)), uni_init(2) + normrnd(0,1.4^2)];
        par_uni = fminsearch(@(par) negL(par, uni_exceed, u_uni), uni_init_temp);
    end
    while par_exp == exp_init_temp
        exp_init_temp = [max(0.1, exp_init(1) + normrnd(0, 1.4^2)), exp_init(2) + normrnd(0,1.4^2)];
        par_exp = fminsearch(@(par) negL(par, exp_exceed, u_exp), exp_init_temp);
    end

    pu = length(uni_exceed)/n;
    p0_uni = pu*(max(0,1 + par_uni(2)*(u0 - u_uni)/par_uni(1)) )^(-1/par_uni(2));

    pu = length(exp_exceed)/n;
    p0_exp = pu*(max(0,1 + par_exp(2)*(u0_trans - u_exp)/par_exp(1)) )^(-1/par_exp(2)); 
    
    % save estimated probabilities
    psave(:,k) = [p0_uni;p0_exp];

end

%% Plots of estimates along with 95% confidence intervals

m_uni = mean(psave(1,:));
m_exp = mean(psave(2,:));
sdev_uni = std(psave(1,:))/sqrt(N);
sdev_exp = std(psave(2,:))/sqrt(N);
ci_uni = m_uni + [-1 1]*1.96*sdev_uni
ci_exp = m_exp + [-1 1]*1.96*sdev_exp
clf
subplot(211)
plot(psave(1,:),'.')
hold on; plot(m_uni*ones(1,N)); plot(ci_uni(1)*ones(1,N)); plot(ci_uni(2)*ones(1,N)); hold off
title('estimate using non-transformed data')
subplot(212)
plot(psave(2,:),'.')
hold on; plot(m_exp*ones(1,N)); plot(ci_exp(1)*ones(1,N)); plot(ci_exp(2)*ones(1,N)); hold off
title('estimate using transformed data')

