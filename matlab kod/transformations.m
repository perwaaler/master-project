%% ANALYSIS ON TRANSFORMED DATA!!!
clf;plot(min(DAFEA'),'.'); ylim([0,20])
%% reciprical transformation analysis
NN = length(DAFEA(:,1));
clf
delta = 0.1
trans_DAFEA = 1./(delta+DAFEA).^5;
plot(log(max(trans_DAFEA')),'.'); hold on
plot(ones(1,NN)*log(.001)); plot(ones(1,NN)*log(.007))
%%
delta = 0.1
k=5;
trans_DAFEA = 1./(delta+DAFEA).^k;
p_EA = (sum(enc_type==-1)+sum(enc_type==-2) + sum(enc_type==2))/N;
m=20;
init = [2 .8];         % initial guess
U = linspace(3e-7 ,3.3546e-04,m);
parameters = zeros(2,m);
p_nea = zeros(1,m);
for k=1:m
    data = trans_DAFEA(:);
    data = data(find(data>U(k)));
    negL = @(par) -sum( log(gppdf(data,par(2),par(1),U(k))) );
    param = fminsearch(negL,init);
    while param == init                                                    % in case initial guess is bad
        init = [max(0.1,init(1) + normrnd(0,1.4^2)), init(2) + normrnd(0,1.4^2)]
        param = fminsearch(negL,init)
    end
    parameters(:,k) = param
    p_u = sum(sum((trans_DAFEA)>U(k)))/( length(DAFEA(1,:))*length(DAFEA(:,1)) )
    ue = U(k) - param(1)/param(2)
    p_nea(k) = p_u*(max(0,1 + param(2)*(1/delta^k - U(k))/param(1)) )^(-1/param(2)) * p_EA
end
log(p_nea)
clf; 
subplot(211)
plot(U,parameters(2,:)) 
subplot(212)
plot(U,p_nea)
%% exponential tranformation
NN = length(DAFEA(:,1));
delta=0.0;
s = .6; 
trans_DAFEA = exp(-s*(DAFEA + delta));
clf; plot(max(trans_DAFEA'),'.'); hold on; plot(ones(1,NN)*0.03); plot(ones(1,NN)*0.25); hold off
%%
delta=-0.0
s = 0.6;
trans_DAFEA = exp(-s*(DAFEA(:) + delta));
p_EA = (sum(enc_type==-1)+sum(enc_type==-2) + sum(enc_type==2))/N;
m=25;
init = [2 .8];         % initial parameter guess
U = linspace(.03 , .25, m);
parameters = zeros(2, m);
p_nea = zeros(1, m);
for k=1:m
    data = trans_DAFEA(:);
    data = data(find(data>U(k)));
    negL = @(par) -sum( log(gppdf(data, par(2), par(1), U(k))) );
    param = fminsearch(negL,init);
    while param == init                                                    % in case initial guess is bad
        init = [max(0.1,init(1) + normrnd(0,1.4^2)), init(2) + normrnd(0,1.4^2)];
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
plot(U,p_nea,'.')
%6e-3 1.5e-5 5e-4 1e-4 .5e-3 1e-3 0.4e-3 .5e-3 .5e-3%%% 4e-4 2e-4 2e-4
%2e-4 6e-5 1.6e-4 3e-4 3e-4 2e-5 8e-4 3e-5  7e-4
%% plotting mah results
clf
plot([9.00E-04	0	4.30E-05	0	4.00E-04	0	0	0	0])
hold on     
plot([1.00E-03	2.10E-03	1.00E-06	0	0	1.50E-04	0	2.40E-03	0])
plot([7.00E-04	8.00E-04	2.00E-04	5.00E-07	1.00E-04	1.00E-04	1.00E-04	2.00E-04	3.00E-05])
plot([6.00E-03	1.50E-05	5.00E-04	1.00E-04	5.00E-04	1.00E-03	4.00E-04	5.00E-04	5.00E-04])







