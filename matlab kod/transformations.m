%% ANALYSIS ON TRANSFORMED DATA!!!
Nenc = length(DAFEA(:,1))
u_lower = 5
u_upper = 10
clf;plot(min(DAFEA'),'.'); hold on; plot(ones(1,Nenc)*u_lower); plot(ones(1,Nenc)*u_upper)
ylim([0,30])
%% reciprical transformation analysis
p=4
delta = 1
%trans = @(x)1./(delta + x).^p
trans = @(x)exp(-1*(x - 0))
%trans = @(x) -x
Nenc = length(DAFEA(:,1));
clf


u_lower_trans = trans(u_upper)
u_upper_trans = trans(u_lower)

trans_DAFEA = trans(DAFEA);
clf; plot(max(trans_DAFEA'),'.'); hold on
plot(ones(1,Nenc)*u_lower_trans); plot(ones(1,Nenc)*u_upper_trans)
%%

Nenc = length(DAFEA(:,1));         % number of encounters for which there was interaction
compute_ci = 1;                    % set equal to one if confidence intervals for xi are desired
Nbs = 100;                         % number of bootstrapped samples to compute standard error
trans_DAFEA = trans(DAFEA);
p_EA = (sum(enc_type==-1)+sum(enc_type==-2) + sum(enc_type==2))/N; % probability for encounter to be interactive
m = 10;                                                    % number of thresholds used for estimation
init = [1 .8];                                             % initial guess
U = linspace(u_lower_trans, u_upper_trans,m);                                 % vector containing thresholds
ci_xi_u = zeros(2,m);                                      % collects 95% ci's for xi
parameters = zeros(2,m);
p_nea = zeros(1,m);                                        % collects estimated collision probability for each threshold

for k=1:m

    data = trans_DAFEA(:);
    data = data(find(data>U(k)));
    negL = @(par) -sum( log(gppdf(data,par(2),par(1),U(k))) );
    param = fminsearch(negL,init);
    while param == init                                                    % in case initial guess is bad
        init = [max(0.1, init(1) + normrnd(0,1.4^2)), init(2) + normrnd(0,1.4^2)]
        param = fminsearch(negL,init)
    end
    parameters(:,k) = param;
    p_u = sum(sum((trans_DAFEA)>U(k)))/( length(DAFEA(1,:))*length(DAFEA(:,1)) );
    p_nea(k) = p_u*(max(0,1 + param(2)*(trans(0) - U(k))/param(1)) )^(-1/param(2)) * p_EA


%%%%%%%%%%%%%%%%%%%%%%%%%% bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if compute_ci == 1
        for j=1:Nbs
            resampling = randsample(Nenc,Nenc,true);              % indeces used to bootstrap
            trans_DAFEA_bs = trans_DAFEA(resampling,:);           % resample
            data = trans_DAFEA_bs(:);
            data = data(find(data>U(k)));       % get exceedences
            negL = @(par) -sum( log(gppdf(data,par(2), par(1), U(k))) ); % negative log likelihood function
            param_bs = fminsearch(negL,param);                           % estimate parameters
            if param_bs == param                                         % in case of stuck
            stuck = 1;
                for t = 1:200
                    t
                    init_temp = [max(0.1,param(1) + normrnd(0.01, 2.0^2) ),param(2) + normrnd(0, 2.0^2)];
                    if negL(init_temp) ~= Inf
                       param_bs = fminsearch(negL,init_temp);
                       if param_bs ~= init_temp
                           break
                       end
                    end

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
        xi_sample(find(xi_sample == NaN)) = [];
        xi_mean = mean(xi_sample);
        se_xi = sqrt(sum((xi_sample - xi_mean).^2)/length(xi_sample));
        ci_xi = param(2) + [-1 1]*1.96*se_xi; ci_xi = sort(ci_xi);
        ci_xi_u(:,k) = ci_xi';
    end

end

clf;
subplot(211)
plot(U,parameters(2,:))
hold on
if compute_ci == 1
    plot(U,ci_xi_u)
end
subplot(212)
plot(U,p_nea)


































%%%%%%%%%%%%%% code from the olden days %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% exponential tranformation
NN = length(DAFEA(:,1));
delta = -0.01;
u_lower_trans = exp(-s*(u_upper)).^p
u_upper_trans = 1/(delta + u_lower).^p
s = .6;
trans_DAFEA = exp(-s*(DAFEA + delta));
clf; plot(mean(trans_DAFEA'),'.'); hold on; plot(ones(1,NN)*0.03); plot(ones(1,NN)*0.2); hold off
%% Estimating using
compute_ci = 1;                                            % set equal to one if confidence intervals for xi are desired
Nbs = 100;                                                 % number of bootstrapped samples to compute standard error
delta= -0.01
s = 0.4;
compute_ci = 1;
trans_DAFEA = exp(-s*(DAFEA + delta));
p_EA = (sum(enc_type==-1)+sum(enc_type==-2) + sum(enc_type==2))/N;
m=10;
ci_xi_u = zeros(2,m);                                      % collects 95% ci's for xi
init = [2 .8];         % initial parameter guess
U = linspace(.01 , .2, m);
parameters = zeros(2, m);
p_nea = zeros(1, m);

for k=1:m
    data = trans_DAFEA(:);
    data = data(find(data>U(k)));
    negL = @(par) -sum( log(gppdf(data, par(2), par(1), U(k))) );
    param = fminsearch(negL,init);
    init_temp = init;
    while param == init_temp                                                    % in case initial guess is bad
        init_temp = [max(0.1,init(1) + normrnd(0,1.4^2)), init(2) + normrnd(0,1.4^2)];
        param = fminsearch(negL,init)
    end
    parameters(:,k) = param
    p_u = sum(sum((trans_DAFEA)>U(k)))/( length(DAFEA(1,:))*length(DAFEA(:,1)) )
    ue = U(k) - param(1)/param(2)
    p_nea(k) = p_u*(max(0,1 + param(2)*(exp(-delta*s) - U(k))/param(1)) )^(-1/param(2)) * p_EA
        %%% bootstrapping %%%%
    if compute_ci == 1
        for j=1:Nbs
            resampling = randsample(Nenc,Nenc,true);  % indeces used to bootstrap
            trans_DAFEA_bs = trans_DAFEA(resampling,:);           % resample
            data = trans_DAFEA_bs(:);
            data = data(find(data>U(k)));       % get exceedences
            negL = @(par) -sum( log(gppdf(data,par(2),par(1),U(k))) ); % negative log likelihood function
            param_bs = fminsearch(negL,param);                          % estimate parameters
            if param_bs == param                       % in case of stuck
            stuck = 1;
                for t = 1:200
                    t
                    init_temp = [max(0.1,param(1) + normrnd(0.01, 2.0^2) ),param(2) + normrnd(0, 2.0^2)];
                    if negL(init_temp) ~= Inf
                       param_bs = fminsearch(negL,init_temp);
                       if param_bs ~= init_temp
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
        xi_sample(find(xi_sample == NaN)) = [];
        xi_mean = mean(xi_sample);
        se_xi = sqrt(sum((xi_sample - xi_mean).^2)/length(xi_sample));
        ci_xi = param(2) + [-1 1]*1.96*se_xi; ci_xi = sort(ci_xi);
        ci_xi_u(:,k) = ci_xi';
    end
end

clf;
subplot(211)
plot(U,parameters(2,:))
hold on
if compute_ci == 1
    plot(U,ci_xi_u)
end
subplot(212)
plot(U,p_nea)














%% exponential tranformation (old code)
NN = length(DAFEA(:,1));
delta=-0.01;
s = .4;
trans_DAFEA = exp(-s*(DAFEA + delta));
clf; plot(mean(trans_DAFEA'),'.'); hold on; plot(ones(1,NN)*0.03); plot(ones(1,NN)*0.1); hold off
%%
delta= -0.01
s = 0.4;
trans_DAFEA = exp(-s*(DAFEA(:) + delta));
p_EA = (sum(enc_type==-1)+sum(enc_type==-2) + sum(enc_type==2))/N;              % probability of evasive action for one encounter
m=5;
init = [2 .8];                     % initial parameter guess
U = linspace(.02 , .1, m);
parameters = zeros(2, m);
p_nea = zeros(1, m);
for k=1:m
    data = trans_DAFEA(:);
    data = data(find(data>U(k)));
    negL = @(par) -sum( log(gppdf(data, par(2), par(1), U(k))) );
    param = fminsearch(negL,init);
    init_temp = init;
    while param == init_temp                                                    % in case initial guess is bad
        init_temp = [max(0.1,init(1) + normrnd(0,1.4^2)), init(2) + normrnd(0,1.4^2)];
        param = fminsearch(negL,init)
    end
    parameters(:,k) = param
    p_u = sum(sum((trans_DAFEA)>U(k)))/( length(DAFEA(1,:))*length(DAFEA(:,1)) )
    ue = U(k) - param(1)/param(2)
    p_nea(k) = p_u*(max(0,1 + param(2)*(exp(-delta*s) - U(k))/param(1)) )^(-1/param(2)) * p_EA
end

clf;
subplot(211)
plot(U,parameters(2,:))
subplot(212)
plot(U,p_nea)
