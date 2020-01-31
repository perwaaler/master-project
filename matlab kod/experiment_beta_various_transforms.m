% The purpose of this code is to show an example of a situation in which
% applying some kind of inversion to the data improves the quality of the
% estimate obtained when using POT method. The data is sampled from a beta
% distribution picked specifically so that the situation mimics that of the
% traffic-security situation; i.e. we are looking to estimate the
% probability of the data falling below a certain threshold, and the
% probability of that happening is extremely low, resulting in the probability
% of dipping below such a threshold often being estimated to be zero.


N = 2000;           % number of iterations
n = 2000;           % number of samples generated each iteration

% parameter values of the beta distribution
a = 6;
b = 5.5;

% initial values when minimizing negtive log-likelihood fcn.
beta_init = [0.1 -0.1];
inv_init = [1 1];                                                                                                 % clf;plot(linspace(0,1,1000), betapdf(linspace(0,1,1000),a,b) )

u0 = 0.08;
p = betacdf(u0,a,b);  % true value of probability we are trying to estimate = 6.0178e-05
psave = zeros(2,N);     % saves estimates for each iteration
par_beta_save = zeros(2,N);
par_inv_save = zeros(2,N);
ue_beta = zeros(1,N);   % saves estimates of lower endpoints of beta distribution
trans = @(x) 1./x.^4    % specifies the type of transformation to use
min_beta = zeros(1,N)   % collects the smallest value from each beta-sample. Used to compare against estimated lower endpoint

u0 = -0.08;             % we are estimating P(X>u0)
u0_inv = trans(-u0);    % corresponding threshold for transformed data
u_beta = -0.4;          % threshold used when applying POT method
u_inv = trans(-u_beta); % corresponding threshold for transformed data
qq_plot = 0;            % set equal to one if qq-plots are desired

for k=1:N
    % generate samples
    beta_sample = betarnd(a*ones(1,n), b*ones(1,n));
    min_beta(k) = min(beta_sample);
    inv_sample = trans(beta_sample);  % inverted data sample
    beta_sample = -beta_sample;       % negate data so that we can apply POT method
    
    % find exceedences
    beta_exceed = beta_sample(find(beta_sample > u_beta));
    inv_exceed = inv_sample(find(inv_sample > u_inv));

    % estimate parameters
    negL = @(par, data, u) -sum( log(gppdf(data,par(2),par(1),u)) );
    beta_init_temp = beta_init;
    inv_init_temp = inv_init;
    par_beta = fminsearch(@(par) negL(par, beta_exceed, u_beta), beta_init);
    par_inv = fminsearch(@(par) negL(par, inv_exceed, u_inv), inv_init);
    
    % in case minimizer gets stuck on  initial guess
    while par_beta == beta_init_temp
        beta_init_temp = [max(0.1, beta_init(1) + normrnd(0,0.8) ), beta_init(2) + norm(0,1)]; 
        par_beta = fminsearch(@(par) negL(par, beta_exceed, u_beta), beta_init_temp);
    end
    while par_inv == inv_init_temp
        inv_init_temp = [max(0.1, inv_init(1) + normrnd(0,0.8) ), inv_init(2) + norm(0,1)]; 
        par_inv = fminsearch(@(par) negL(par, inv_exceed, u_inv), inv_init_temp);
    end
    
    % estimate P(X>x0)
    pu = length(beta_exceed)/n;
    p0_beta = pu*(max(0,1 + par_beta(2)*(u0 - u_beta)/par_beta(1)) )^(-1/par_beta(2));

    pu = length(inv_exceed)/n;
    p0_inv = pu*(max(0,1 + par_inv(2)*(u0_inv - u_inv)/par_inv(1)) )^(-1/par_inv(2));
    
    % save estimates
    psave(:,k) = [p0_beta;p0_inv];
    par_beta_save(:,k) = par_beta';
    par_inv_save(:,k) = par_inv';
    
    % endpoints
    ue_beta(k) = u_beta - par_beta(1)/par_beta(2);
    
    % qq plots
     if qq_plot == 1
        F_inv = @(y,scale,shape) scale/shape*((1-y).^-shape - 1);
        N_inv = length(inv_exceed);
        N_beta =  length(beta_exceed);
        inv_excess = inv_exceed - u_inv;
        beta_excess = beta_exceed - u_beta;
        clf 
        subplot(211)
            plot(sort(beta_excess), F_inv( [1:N_beta] /(N_beta+1), par_beta(1),par_beta(2)) ,'.');
            hold on
            plot(sort(beta_excess), sort(beta_excess)); hold off 
            title('no transformation')
        subplot(212)
            plot(sort(inv_excess),F_inv( [1:N_inv] /(N_inv+1), par_inv(1),par_inv(2)) ,'.');
            hold on
            plot(sort(inv_excess), sort(inv_excess)); hold off
            title('transformation')
        pause(0.5)
     end
    
     beta_init = par_beta;
     inv_init = par_inv;  
end

%% Plots
m_beta = mean(psave(1,:));
m_inv = mean(psave(2,:));
avg_error_beta = sqrt(sum((psave(1,:)-p).^2)/N)
avg_error_inv = sqrt(sum((psave(2,:)-p).^2)/N)
clf
subplot(211)
    plot(psave(1,:),'.'); hold on 
    plot(ones(1,N)*p, 'k')
    plot(ones(1,N)*m_beta)
    title('estimate using non-transformed data')
    ylim([0,3e-3])
subplot(212)
    plot(psave(2,:),'.'); hold on
    title('estimate using transformed data')
    plot(ones(1,N)*p,'k'); plot(ones(1,N)*m_inv)
    ylim([0,3e-3])
%% plot of estimated lower endpoints of beta distribution
clf
subplot(211)
plot(linspace(0,1,1000), betapdf(linspace(0,1,1000),a,b) )
title('pdf of beta distribution that generated the data')
subplot(212)
plot(-ue_beta,'.'); hold on
plot(min_beta,'.'); hold off
title('estimates(blue) of lower endpoints of beta distribution')
ylim([-1,1])
