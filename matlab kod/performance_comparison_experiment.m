
%% Estimation of P(collision)
% this code is used to take simulated data and perform POT method to
% estimate P(collision). Data_matrix  can be either danger_FEA, DAFEA, X,
% or danger_max_EA. data_type=1 --> danger_FEA or DAFEA, data_type=2 --> X,
% data_type=3 --> danger_max_EA


data_type = 2
L = size(all_data,1);
m = 10;                                                                % number of thresholds used for estimation
ue_save_matrix = zeros(L,m);
param_save_matrix = zeros(L,m);
pc_save_matrix = zeros(L,m);
p_c_save_matrix = zeros(L,m);
p_interactive_matrix = all_data{:,7}
p_ea_matrix = all_data{:,8}


% transformation parameters
select_trans = 1;
p_ex = 0.6;
p_inv = 2;
d_inv = 3;

% plotting options
compute_ci = 0;                                       % set equal to one if confidence intervals for xi are desired
qqplot = 0;
save_plot = 0;

% pausing options
pause_trans = 0;
qq_pause = 0.0;
stability_pause = 0.0

for l=1:L
% generate sample of encounters
l
data_matrix = all_data{l,data_type}; % select data. row i should correspond to encounter i, and column j to j'th simulated ttc value (in case of stochastic ttc)


% find encounters with finite ttc values
min_data = data_matrix(find(min(data_matrix,[],2)<Inf),:);
min_data = min(min_data,[],2);

% find minimum and maximum thresholds based on amount of data to be used
u_minmax = find_threshold(min_data, 0.2, 0.5);
u_l = u_minmax(1);
u_u = u_minmax(2);

clf;plot(min_data,'.'); hold on; plot(ones(1,length(min_data))*u_l); plot(ones(1,length(min_data))*u_u)
ylim([-1,30])

%%% transforming data and plotting transformed data and thresholds
% transforms
transinv = @(x)1./(d_inv + x).^p_inv;
transex = @(x)exp(-p_ex*(x - 2));
transneg = @(x) -x;

% choose transform to use
if select_trans == 1
    trans = @(x) transneg(x);
elseif select_trans == 2
    trans = @(x) transex(x);
elseif select_trans == 3
    trans = @(x) transinv(x);
end

Nenc = length(data_matrix(:,1));

% transform thresholds
u_min_max = sort(trans([u_l,u_u]));
u_l_trans = u_min_max(1); 
u_u_trans = u_min_max(2);
m = 10;
U = sort(trans(linspace(u_l, u_u,m)));

trans_data = trans(data_matrix);

clf; plot(trans(min_data),'.'); hold on
plot(1:length(min_data),U'*ones(1,length(min_data)), 'k');plot(ones(1,length(min_data))*trans(0));ylim([trans(30),trans(0)*1.1])
title('transformed data')
pause(pause_trans)

%%% ensure good initial value!
trans_data = trans_data(:);
exceed = trans_data(find(trans_data > u_l_trans));
negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
init = fminsearch(@(par) negL(par, exceed, u_l_trans), [3 0.3]);
pause(1)
%%% Fit models for range of thresholds and show diagnostic plots
shake_guess = 0.1;                 % variance of noise that gets added to initial guess when stuck
Nenc = length(data_matrix(:,1));   % number of encounters
Nexp = length(data_matrix(1,:));   % number of values per row. If surrogate-measure is deterministic, then Nexp=1.
Nbs = 100;                         % number of bootstrapped samples to compute standard error
trans_data = trans(data_matrix);   
U = sort(trans(linspace(u_l, u_u,m)));                                 % vector containing thresholds'

% collects 95% ci's for xi
ci_sigma_u = zeros(2,m)*nan;
ci_xi_u = zeros(2,m)*nan;                                      % collects 95% ci's for xi
ci_p_nea_u = zeros(2,m)*nan;

param_save = zeros(2,m);
pc = zeros(1,m)*nan;                                    % collects estimated collision probability for each threshold
ue_save = zeros(1,m)*nan;                                  % collects estimated upper endpoint
max_data = max(max(trans_data));                           % largest observed value
negL = @(par, exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );  %negative log likelihood fcn.
logit = 1;                                               % set to plot logarithm of p_nea when magnitude of p_nea varies alot


% limits for plots of empirical and model distribution functions
xplot_lower = 0;
xplot_upper = 1;
n_eval_cdf = 1000;    % number of points where cdf gets evaluated

for k=1:m
    k;
    data = trans_data(:);
    exceed = data(find(data>U(k)));
    param = fminsearch(@(par) negL(par,exceed,U(k)),init);
    init_temp = init;
    while_counter = 0;
    while param == init_temp                                                    % in case initial guess is bad
        while_counter = while_counter + 1;
        init_temp = [max(0.1, init(1) + normrnd(0,shake_guess^2)), init(2) + normrnd(0,shake_guess^2)];
        param = fminsearch(@(par) negL(par,exceed,U(k)),init_temp);
        if while_counter == 10
            break
            'pick better initial value'
        end
    end
    param_save(:,k) = param;
    p_u = length(exceed)/length(data);
    pc(k) = p_u*(1 - gpcdf(trans(0), param(2), param(1),U(k)) );
    ue = U(k) - param(1)/param(2);
    if param(2)<0; ue_save(k) = ue; end 
    
    if qqplot == 1;
        
        %%% generate sample of stochastic ttc, where one sample is drawn from each encounter.                   
        col_ind = randsample(Nexp,Nenc, true)';
        ind = sub2ind(size(trans_data), 1:Nenc, col_ind );
        ttc_sample = trans_data(ind); 
        exceed = ttc_sample(find(ttc_sample>U(k)));
        excess = exceed - U(k);
        
        %%% compute empirical d.f. for excesses
        exceed_data = trans_data(    find(max(trans_data,[],2) > U(k)) ,:   );   % keeps only rows with atleast one obs. above threshold
        pu_i = sum(exceed_data>U(k), 2)/ size(trans_data,2);                     % probability to exceed threshold for each encounter
        weights = pu_i/sum(pu_i) ;                                          

        %%% replace obs. below u with nan's
        below_u = find(exceed_data<=U(k));
        exceed_data(below_u) = nan*ones(1,length(below_u));
        exceed_data = exceed_data - U(k); % subtract u to get excessses

        %%% evaluate empirical distribution function
        xplot_lower = 0;
        xplot_upper = max(excess)*3;
        x_eval = linspace(xplot_lower, xplot_upper, n_eval_cdf);
        femp = weights'*F_emp(x_eval, exceed_data);
        
        %%% plot distribution functions and qq-plots
        clf
        subplot(211)
            qq_plot(exceed,param(1),param(2),U(k),k)
            title(sprintf('empirical vs model quantiles, for threshold %d.',k))
        subplot(212)
            plot(x_eval,femp,'.')
            hold on 
            plot(x_eval, gpcdf(x_eval, param(2), param(1), 0))
            line(([0,0]), 1.2,'LineStyle','--');
            line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');
            title(sprintf('empirical d.f. vs model d.f. for threshold %d.',k))
        % save plot
        if save_plot == 1
            savefig(sprintf('goodness_of_fit_test_DAFEA_exp_%d_%d_id4',l,k))
        end
%             histogram(excess,'Normalization','probability'); hold on
%             plot(linspace(0,max(excess)*1.5,100), gppdf(linspace(0,max(excess)*1.5,100),param(2),param(1),0) )
%             hold off
        pause(qq_pause)
    end

    init = param;
    
%%%%%%%%%%%%%%%%%%%%%%%%%% bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if compute_ci == 1
        
        for j=1:Nbs
            
            resampling = randsample(Nenc,Nenc,true);              % indeces used to bootstrap
            trans_DAFEA_bs = trans_data(resampling,:);            % bootstrapped sample
            data = trans_DAFEA_bs(:);
            exceed = data(find(data>U(k)));                         % get exceedences
            param_bs =  fminsearch(@(par) negL(par,exceed,U(k)),init);                           % estimate parameters
            init_temp = init;
            while_counter = 1;
            
            while param_bs == init_temp 
                while_counter = while_counter + 1;
                % in case initial guess is bad
                init_temp = [max(0.1, init(1) + normrnd(0,shake_guess^2)), init(2) + normrnd(0,shake_guess^2)];
                param_bs = fminsearch(@(par) negL(par,exceed,U(k)),init_temp);
                if while_counter==30; param_bs = [nan,nan]; break; end
            end
            
            p_u = length(exceed)/length(data);
            sigma_sample(j) = param_bs(1);
            xi_sample(j) = param_bs(2);
            p_nea_sample(j) = p_u*(1 - gpcdf(trans(0), param_bs(2), param_bs(1),U(k)) );
            
        end
        sigma_sample(find(isnan(sigma_sample))) = [];
        xi_sample(find(isnan(xi_sample))) = [];
        p_nea_sample(find(isnan(p_nea_sample))) = [];
        
        sigma_mean = mean(sigma_sample);
        xi_mean = mean(xi_sample);
        p_nea_mean = mean(p_nea_sample);
        
        se_sigma = sqrt(sum((sigma_sample - sigma_mean).^2)/length(sigma_sample));
        se_xi = sqrt(sum((xi_sample - xi_mean).^2)/length(xi_sample));
        se_p_nea = sqrt(sum((p_nea_sample - p_nea_mean).^2)/length(p_nea_sample));
        
        ci_sigma = param(1) + [-1 1]*1.96*se_sigma; ci_sigma = sort(ci_sigma);
        ci_xi = param(2) + [-1 1]*1.96*se_xi; ci_xi = sort(ci_xi);
        ci_p_nea = pc(k) + [-1 1]*1.96*se_p_nea; ci_p_nea = sort(ci_p_nea);
        
        ci_sigma_u(:,k) = ci_sigma';
        ci_xi_u(:,k) = ci_xi';
        ci_p_nea_u(:,k) = ci_p_nea';
        
    end
    
end

clf;
subplot(221)
plot(U,param_save(1,:))
title('sigma_{est}')
hold on
if compute_ci == 1
    plot(U,ci_sigma_u)
end
subplot(222)
plot(U, param_save(2,:))
title('xi_{est}')
hold on
if compute_ci == 1
    plot(U,ci_xi_u)
end

subplot(223)
if min(param_save(2,:))<0
    plot(U,ue_save); hold on
    plot(U,ones(1,m)*max_data)
    plot(U,ones(1,m)*trans(0))
    title('upper endpoint estimates')
end
subplot(224)    
if logit==1
    plot(U,log10(pc))
    title('log(p_{est})')
    if compute_ci == 1
        hold on
        %plot(U,ci_p_nea_u)
    end
else
    plot(U,pc)
    title('p_{est}')
    if compute_ci == 1
        hold on
        %plot(U,ci_p_nea_u)
    end
end

if save_plot == 1
    savefig(sprintf('range_of_threshold_plots_DAFEA_exp_%d_id4',l))
end


param_save_matrix(l,:) = param_save(1,:) + 1i*param_save(2,:);
ue_save_matrix(l,:) = ue_save;
pc_save_matrix(l,:) = pc;
pause(stability_pause)


%%% estimating p(collision type i)
p_interactive = all_data{l,7}; % probability of encounter being interactive
p_ea = all_data{l,8};
% this estimator gives P(C1,NEA) when p_nea is based on danger at first
% interactive action, DAFEA
if data_type == 1
    p_c = p_interactive*pc
end

% this estimaor gives P(C, NEA) when p_nea is based on X
if data_type == 2
    p_c = pc
end

% this estimaor gives P(C, EA) when p_nea is based on danger_FEA
if data_type == 3
    p_c = pc*all_data{l,7}
end

% this estimaor gives P(C, EA) when p_nea is based on danger_max_EA
if data_type == 4
    p_c = pc*all_data{l,7}
end


p_c_save_matrix(l,:) = p_c;

end
%% save information
hit_rate = sum(pc_save_matrix>0,1)/500*100;
%%
save('hit_rate_neg_X_id5','hit_rate')
save('pc_neg_X_id5', 'pc_save_matrix')
save('p_neg_X_id5', 'p_c_save_matrix')
save('param_neg_X_id5', 'param_save_matrix')
save('ue_neg_X_id5','ue_save_matrix')

