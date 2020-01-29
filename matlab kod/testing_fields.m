for k=1:length(DAFEA(1,:))
    plot(DAFEA(k,:),'.')
    pause(0.3)
end

%% cleaning away all rows wich correspond to clean crashes:
DAFEA(find(enc_type==1),:) = [];
%%
datamatrix = DAFEA;
m = 25;
N_ea = length(DAFEA(:,1)); % number of encounters for which there was evasive action
U = linspace(-5,-1.5,m);
parameters = zeros(2,m);
p_nea = [1,m];
for k = 1:m;
f = @(param)negloglik(param,U(k),N_ea,datamatrix);
param = fminsearch(f,[2 1])
parameters(:,k) = param;
p_u = sum(sum(DAFEA'<-U(k)))/(NTTC*N_ea);
p_nea(k) = p_u*(max(0,1 + param(2)*(0-U(k))/param(1)) )^(-1/param(2))
end

plot(U,parameters(1,:) - U.*parameters(2,:)); ylim([-1 1])
plot(U,p_nea); ylim([-1*10^-3,0.3*10^-2])
% pest  0.00029345 0.000329 0.00036455 
% 0.8287, 0.0012
% matlab: 0.0019847 R: 0
%%
datamatrix2 = danger_FEA';
m = 1;
N_ea2 = length(datamatrix2);
U2 = linspace(-2,-1,m);
parameters2 = zeros(2,m);
p_nea2 = [1,m];
for k = 1:m
    negL = @(par) -sum( log(gppdf(negdata,par(2),par(1),U2(k))) );
    param2 = fminsearch(negL,[0.5684280 -0.6288933])
    parameters2(:,k) = param2
    p_u2 = sum(danger_FEA'<-U2(k))/(N_ea2);
    p_nea2(k) = p_u2*(max(0,1 + param2(2)*(0-U2(k))/param2(1)) )^(-1/param2(2));
end

plot(U2,parameters2(1,:) - U2.*parameters2(2,:)); ylim([-1 1])
plot(U2,p_nea2); ylim([-1*10^-3,0.3*10^-2])



negL = @(par) -sum( log(gppdf(negdata,par(2),par(1),U2(k))) );
param2 = fminsearch(negL,[0.5684280 -0.6288933])
