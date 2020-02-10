% code for computing and plotting empirical vs model d.f.

k=9
xplot_lower = 0;
xplot_upper = 2.5;
u = U(k);                                                             % set desired threshold
exceed_data = trans_data(    find(max(trans_data,[],2) > u) ,:   ); % matrix with data only from encounters with atleast one exceedence of threshold
pu_i = sum(exceed_data>u, 2)/ size(trans_data,2);                   % probability to exceed threshold for each encounter
weights = pu_i/sum(pu_i) ;                                           %

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

clf
plot(linspace(xplot_lower,xplot_upper,1000),weights'*F_emp(linspace(xplot_lower,xplot_upper,1000),exceed_data),'.')
hold on
plot(linspace(xplot_lower,xplot_upper,1000), gpcdf(linspace(xplot_lower,xplot_upper,1000), param(2),param(1),0))
line(trans([0,0]), get(gca, 'ylim'),'LineStyle','--');
line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');