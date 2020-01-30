u = -3.0
N_ea = length(DAFEA(:,1));             % row number of DAFEA
min_ttc = zeros(1,N_ea);
mean_ttc = zeros(1,N_ea);
index_mat = zeros(size(DAFEA));
pot = 0;
index_EA = [];
for n=1:N_ea
    ttcvec = DAFEA(n,:);
    index_mat(n,:) = ttcvec < 1000;
    ttcvec = ttcvec(find(ttcvec < 1000));
    min_ttc(n) = min(DAFEA(n,:));
    mean_ttc(n) = mean(ttcvec);
    ttcvec = [ttcvec, zeros(1, NTTC-length(ttcvec))];
    pot = pot + sum(DAFEA(n,:)<-u);
end
pu = pot/(N*N_ea)

datamatrix=DAFEA;
% U = linspace(-60,-30,15);
% parameters = zeros(2,5);
% for k = 1:length(U);
% f = @(param)negloglik(param,U(k),N,datamatrix);
% 
% %options = optimset('PlotFcns',@optimplotfval)
% param = fminsearch(f,[2 1])%,options)
% parameters(:,k) = param;
% end
%u - 8.6346 /-0.3862
%param = [9.2780   -0.3916];
%u=-35;
f = @(param)negloglik(param,u,N_ea,datamatrix);
options = optimset('PlotFcns',@optimplotfval)
param = fminsearch(f,[2 1])%, options);
sigma = param(1)
xi = param(2)

pu*(max(0,1 + xi*(0-u)/sigma)  )^(-1/xi)


% 22:0.0052 21:0.0055 20:0.0046 19:0.0047 18:0.0050
%results for progress-presentation: 0.0038 0.0027 0.0016 0.0053 0.0027 0.0023 
% N=400;NTTC=3000; 0.0022 0.0024 0.0045 0.0051

% matlab:0.00026539 R:0.003792645

