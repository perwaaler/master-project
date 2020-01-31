%% optimization for particular threshold u for DAFEA
Nenc = length(DAFEA(:,1));
clf;plot(min(DAFEA'),'.');ylim([-1,10]); hold on; plot(ones(1,Nenc)*5)
%%
Nenc = length(DAFEA(:,1));                                 % number of encounters where there was evasive action
Nbs = 100;                                                 % number of bootstrapped samples
init = [0.3 -.5]                                           % initial parameter guess
negdata = -DAFEA(:);                                       % negate data
m = 10                                                     % number of thresholds
U = linspace(-6,-4,m);                                     % range of thresholds
ci_xi_u = zeros(2,m);                                      % collects 95% ci's for xi
p_nea = zeros(1,m);                                        % vector that collects the estimated probability of pure collision for each threshold
parameters = zeros(2,m);                                   % saves parameter estimates for each threshold
compute_ci = 1;
for k=1:m
    u = U(k);                                                  % initial parameter guess
    negdata = -DAFEA(:);
    negdata = negdata(find(negdata>u));                        % selects exceedences only
    negL = @(par) -sum( log(gppdf(negdata,par(2),par(1),u)) ); % define negative log-likelihood function
    param = fminsearch(negL,init);                             % obtain parameter estimate

    xi_sample = zeros(1,Nbs);
    init_temp = init;

    while param == init_temp                                   % in case initial guess was bad
        init_temp = [max(0.1,init(1) + normrnd(0.1,1^2)), init(2) + normrnd(0,0.7^2)];
        param = fminsearch(negL,init_temp);
    end
    parameters(:,k) = param
    init = param;
    p_u = sum(sum((DAFEA)<-u))/( length(DAFEA(1,:))*length(DAFEA(:,1)) ) % probability to exceed threshold
    %if param(2)<0; ue = u - param(1)/param(2); else ue =Inf; end         % upper endpoint estimate
    p_nea(k) = p_u*(max(0,1 + param(2)*(0 - u)/param(1)) )^(-1/param(2))     % estimate probability of pure collision
    if compute_ci == 1
    %%% bootstrapping %%%%
    for j=1:Nbs
        resampling = randsample(Nenc,Nenc,true);  % indeces used to bootstrap
        DAFEA_bs = DAFEA(resampling,:);           % resample
        negdata = -DAFEA_bs(:);                   % negate data
        negdata = negdata(find(negdata>u));       % get exceedences
        negL = @(par) -sum( log(gppdf(negdata,par(2),par(1),u)) ); % negative log likelihood function
        param_bs = fminsearch(negL,param);                          % estimate parameters
        if param_bs == param                       % in case of stuck                                
        stuck = 1;
            for t = 1:200
                t
                init_temp = [max(0.1,param(1) + normrnd(0.01, 2.0^2) ),param(2) + normrnd(0, 2.0^2)];
                if negL(init_temp) ~= Inf
                   param_bs = fminsearch(negL,init_temp); 
                   if param_bs ~= init_temp
                       stuck
                       break
                   end
                end
%                 if t>195
%                     pause(0.5)
%                 end
                

                if t==200                          % give up if after 100 iterations no results have been obtained
                    param_bs = [NaN,NaN];
                    SAMPLE = DAFEA_bs;
                    PARAM = param;
                    u_corrupt = u;
                end
            end
        end
        xi_sample(j) = param_bs(2);
    end
    end
    xi_sample(find(xi_sample == NaN)) = [];
    xi_mean = mean(xi_sample);
    se_xi = sqrt(sum((xi_sample - xi_mean).^2)/length(xi_sample));
    ci_xi = param(2) + sign(param(2))*[-1 1]*1.96*se_xi
    ci_xi_u(:,k) = ci_xi';
end
%% plots
clf; 
subplot(211)
plot(U,parameters(2,:)) 
hold on
if compute_ci == 1
    plot(U,ci_xi_u)
end
subplot(212)
plot(U,p_nea)
%% investigating difficult case
u =  u_corrupt;
negdata = -SAMPLE(:);                     % negate data
negdata = negdata(find(negdata>u));       % get exceedences
negL = @(par) -sum( log(gppdf(negdata,par(2),par(1),u)) ); % negative log likelihood function
%param_bs = fminsearch(negL,[1 3]);
param_temp = [max(0.1,PARAM(1) + normrnd(0.1, 2.0^2)),PARAM(2) + normrnd(0, 2.0^2)]
negL(param_temp)
options = optimset('PlotFcns',@optimplotfval);
param_bs = fminsearch(negL,param_temp, options)
param_bs==param_temp
