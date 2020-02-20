
data = DAFEA;
trans =@(x) -x;
trans_data = trans(data);
u = 7
exceed_data = data( min(data,[],2) < u ,: );   % keeps only rows with atleast one obs. above threshold
n_exceed = size(exceed_data,1);
gampar_save = zeros(2,n_exceed);
point_masses = [];
for i=1:n_exceed
    sample_i = exceed_data(i,:);
    sample_i = sample_i(sample_i<Inf);
    if length(sample_i)==1
        point_masses(end+1) = sample_i;
    end
    gampar_save(:,i) = gamfit(sample_i)'
    clf
    plot(linspace(0,10,1000), gampdf(linspace(0,10,1000), gampar_save(1,i),gampar_save(2,i)) )
    hold on;
    histogram(sample_i)
    pause(0)
end
%%

plot(trans_data(:),'.')

gampar_save(:,find(gampar_save(1,:)==inf))=[];

integrand =@(x,gpdpar) gampdf(-x,gampar_save(1,1),gampar_save(2,1)) .* log(gppdf(x-(-u),gpdpar(2),gpdpar(1),0));
integral_i =@(gpdpar) integral(@(x) gampdf(-x,gampar_save(1,1),gampar_save(2,1)) .* log(gppdf(x-(-u),gpdpar(2),gpdpar(1),0))  , -u, 0 );




gpdpar = [1 0.01]
summer = 0;
for i=1:length(gampar_save)
summer + -integral(@(x) integrand(x, gampar_save(:,i),-u, gpdpar), -u, 0)
end




