%% plot used to explain danger_FEA and danger_max_EA
data_matrix = all_data{8,3}; % select data. row i should correspond to encounter i, and column j to j'th simulated ttc value (in case of stochastic ttc)
danger_max_EA = all_data{8,4};

% find encounters with finite ttc values
min_data = data_matrix(find(min(data_matrix,[],2)<Inf),:);
min_data = min(min_data,[],2);
coll_index = find(min_data<0);

% find minimum and maximum thresholds based on amount of data to be used
u_minmax_FEA = find_threshold(min_data, 0.2, 0.5);
u_minmax_EA = find_threshold(danger_max_EA,0.2,0.5)
u_FEA = u_minmax_FEA(1);
u_EA = u_minmax_EA(1);

clf;
subplot(211)
plot(min_data,'.'); hold on; 
plot(ones(1,length(min_data))*u_FEA,'--');
plot(ones(1,length(min_data))*0,'--')
plot(coll_index,min_data(coll_index),'*','color','b')
ylim([-1,5])
title('separation-distance at moment of first evasive action')
xlabel('encounter')
ylabel('separation-distance')
subplot(212)

plot(danger_max_EA,'.'); hold on
plot(ones(1,length(danger_max_EA))*u_EA,'--');
plot(ones(1,length(danger_max_EA))*0,'--')
ylim([-0.5,3])
title('separation-distance minimum during attempt to avoid collision')
xlabel('encounter')
ylabel('separation-distance')
%% plots used to show the performance 

%% compare p_pure_coll (DAFEA and danger_FEA) from experiment 4 and 5, exp_par = 0.6
% DAFEA, exp0.6
load p_c_ex_DAFEA_6_id4.mat
pc_pure_id4 = p_c_save_matrix;
load p_c_ex_DAFEA_6_id5.mat
pc_pure_id5 = p_c_save_matrix;
success_rate_DAFEA = sum(pc_pure_id5>pc_pure_id4)/500;
fail_rate_DAFEA = sum(pc_pure_id5<pc_pure_id4)/500

% compare p_pure_coll (danger_FEA) from experiment 4 and 5
load p_c_ex_danger_FEA_6_id4.mat
pc_pure_id4 = p_c_save_matrix;
load p_c_ex_danger_FEA_6_id5.mat
pc_pure_id5 = p_c_save_matrix;
success_rate_danger_FEA = sum(pc_pure_id5>pc_pure_id4)/500;
fail_rate_danger_FEA = sum(pc_pure_id5<pc_pure_id4)/500;

% X, exp(-0.6x)
load p_c_ex_X_6_id4.mat
pc_pure_X_id4 = p_c_save_matrix;
load p_c_ex_X_6_id5.mat
pc_pure_X_id5 = p_c_save_matrix;
fail_rate_X = sum(pc_pure_X_id5<pc_pure_X_id4)/500;
success_rate_X = sum(pc_pure_X_id5>pc_pure_X_id4)/500;

% DAFEA negated
load p_c_neg_DAFEA_id4.mat
pc_pure_DAFEA_neg_id4 = p_c_save_matrix;
load p_c_neg_DAFEA_id5.mat
pc_pure_DAFEA_neg_id5 = p_c_save_matrix;
fail_rate_DAFEA_neg = sum(pc_pure_DAFEA_neg_id5<pc_pure_DAFEA_neg_id4)/500;
success_rate_DAFEA_neg = sum(pc_pure_DAFEA_neg_id5>pc_pure_DAFEA_neg_id4)/500;


ci_success_DAFEA = success_rate_DAFEA' +  (sqrt( success_rate_DAFEA.*(1-success_rate_DAFEA)/500 )*1.96)'*[-1 1]
ci_success_danger_FEA = success_rate_danger_FEA' +  (sqrt( success_rate_danger_FEA.*(1-success_rate_danger_FEA)/500 )*1.96)'*[-1 1]
ci_fail_DAFEA = fail_rate_DAFEA' +  (sqrt( fail_rate_DAFEA.*(1-fail_rate_DAFEA)/500 )*1.96)'*[-1 1]
ci_fail_danger_FEA = fail_rate_danger_FEA' +  (sqrt( fail_rate_danger_FEA.*(1-fail_rate_danger_FEA)/500 )*1.96)'*[-1 1]
ci_fail_X = fail_rate_X' + sqrt( fail_rate_X.*(1-fail_rate_X)/500 )'*[-1 1]*1.96
ci_success_X = success_rate_X' + sqrt( success_rate_X.*(1-success_rate_X)/500 )'*[-1 1]*1.96
ci_fail_DAFEA_neg = fail_rate_DAFEA_neg' + sqrt( fail_rate_DAFEA_neg.*(1-fail_rate_DAFEA_neg)/500 )'*[-1 1]*1.96
ci_success_DAFEA_neg = success_rate_DAFEA_neg' + sqrt( success_rate_DAFEA_neg.*(1-success_rate_DAFEA_neg)/500 )'*[-1 1]*1.96


clf
subplot(211)
a = plot(success_rate_DAFEA,'r');hold on
plot(ci_success_DAFEA,'--','color','r')

b = plot(success_rate_danger_FEA,'b')
plot(ci_success_danger_FEA,'--','color','b')

c = plot(success_rate_X,'g')
plot(ci_success_X,'--','color','g')

d = plot(success_rate_DAFEA_neg,'m')
plot(ci_success_DAFEA_neg,'--','color','m')

title('probability to successfully guess the safer traffic location for different methods')
xlabel('threshold number')

subplot(212)
a = plot(fail_rate_DAFEA,'r');hold on
plot(ci_fail_DAFEA,'--','color','r')
b = plot(fail_rate_danger_FEA,'b')
plot(ci_fail_danger_FEA,'--','color','b')
c = plot(fail_rate_X,'g')
plot(ci_fail_X,'--','color','g')
d = plot(fail_rate_DAFEA_neg,'m')
plot(ci_fail_DAFEA_neg,'--','color','m')
legend([a b c d],'T(DAFEA), stoch. TTC','T(DAFEA), sep. dist.', 'T(DAFEA+maxDDEA), stoch. TTC.','DAFEA, stoch. TTC, no transform')
title('probability to incorrectly guess the safer traffic location for different methods')
xlabel('threshold number')
