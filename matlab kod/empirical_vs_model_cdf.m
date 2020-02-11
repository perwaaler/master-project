% code for computing and plotting empirical vs model d.f.

k=7
xplot_lower = 0;
xplot_upper = 2.5;
data_string = trans_data(:);
ordered_excesses = sort(data_string(find(data_string>u))) - u;
u = U(k);                                                             % set desired threshold
exceed_data = trans_data(    find(max(trans_data,[],2) > u) ,:   );   % matrix with data only from encounters with atleast one exceedence of threshold
pu_i = sum(exceed_data>u, 2)/ size(trans_data,2);                     % probability to exceed threshold for each encounter
weights = pu_i/sum(pu_i) ;                                            % weights for each distribution function.

% replace obs. below u with znan's
below_u = find(exceed_data<=u);
exceed_data(below_u) = nan*ones(1,length(below_u));

exceed_data = exceed_data - u;


femp = weights'*F_emp(linspace(0,3,1000), exceed_data);


% estimate model parameters
trans_data_copy = trans_data(:);
exceed = trans_data_copy(find(trans_data_copy(:) > u_l_trans));
negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
param = fminsearch(@(par) negL(par, exceed, u), [param_save(1,k),param_save(2,k)]);

% plot empirical cdf vs model cdf
clf
x = linspace(xplot_lower,xplot_upper,1000)% ordered_excesses;% 
plot(x,weights'*F_emp(x,exceed_data),'.')
hold on
plot(x, gpcdf(x, param(2),param(1),0))
line(trans([0,0]), [0 1.2],'LineStyle','--');
line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');
title('empirical cdf vs model cdf')