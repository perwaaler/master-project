clear all
pause_length = 2^-6;
plotting = 1; % set to one if plots of encounters are wanted

% parameters to determine distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.30*2;
b = 0.2/2;

% parameters for the gamma distribution for the initial steps
a_init = 0.08*10;
b_init = 10^-1;

% variance of the step-size
var_step = 0.0001;

% parameters relating to detection probability, p(d) = exp(-s*(D + d)^p)
A = 1;     % amplitude; reduce to make pd smaller everywhere
d = .60;   % acts as a focal point; pd for D larger than this are suppressed by larger p, while D smaller than this are increased by larger p
s = 0.4;   % decrease to make it harder to detect; lower to make difference smaller for large and small seperation distances D
p = 1.3;   % This parameter allows you to regulate difference in pd for large and small D; if p is small then pd becomes more uniform
% plot(linspace(0,10,100),exp(-s*(linspace(0,10,100)+d).^p))


thetamod = 1.2;                  % determines how quickly the person can change direction
thetavar = 0.1;                  % determins the amount of random variability of the step-angle at each iteration.
r = 0.06;                        % collision radius of each person
xinit = 10;                      % determines how far apart they are at the start
sigma = 0.3;                     % variance of initial starting position
N = 400;                         % number of encounters
NTTC = 3000;                     % Number of TTC to sample at first evasive action for each encounter
danger_FEA = nan*zeros(1,N);            % danger index at time of first evasive action
DAFEA = zeros(N, NTTC);          % Saves danger at first moment of evasive action; each row contains NTTC sampled ttc values for one encounter
danger_max_nodetec = nan*zeros(1,N);    % highest index of danger for non-detection encounters (i.e. type 1 encounters)
danger_nodetec = zeros(1,N);     % highest index of danger for non-detection encounters
enc_type = zeros(1,N);           % vector containing encounter type of each encounter
compute_ttc = 0;                 % tells the algorithm whether or not you want it to estimate ttc distribution for each encounter
min_stepsize = 0.02;             % set to make it impossible for a stepsize to be smaller than this
for i=1:N
i
%%% initiation of encounter
first_detection = 0;     % They have not yet detected eachother
stepsize = min_stepsize*[1 1] + [gamrnd(a_init,b_init),gamrnd(a_init,b_init)]; % determines average speed
A_real = -xinit;
B_real = xinit;
A_im = normrnd(0,sigma);
B_im = normrnd(0,sigma);
A = A_real+1i*A_im;
B = B_real+1i*B_im;
% data collection variables
detection_status = 0; 
encounter_classifier = 1; % 1 --> (no detection, no crash), -1 --> (detection, no crash) 2 --> (no detection, crash), -2 --> (detection, crash)
danger_vec = [];          % collects the danger index associated with each timestep

    while real(A) < xinit && real(B) > -xinit
        D = norm(A-B);
        pd = exp(-s*(D+d)^p); % detection probability
        %%%%%%%%%%%%%%%% detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  (detection_status==1 | rand(1) < pd) && real(A) < real(B)
            % saving information in vectors, and updating status
            if first_detection == 0 & compute_ttc == 1 % a loop that collects information at first moment of evasive action
                for k=1:NTTC
                    ttc = ttc_simulator_momentum_V3(A,B,stepsize,min_stepsize,var_step,r, thetavar);
                    DAFEA(i,k) = ttc;
                end
                danger_FEA(i) = norm(A0-B0) - 2*r;
            end
            first_detection = 1;
            detection_status = 1;
            encounter_classifier = -1;
            danger_index = D - 2*r;
            danger_vec(length(danger_vec) + 1) = danger_index;

            % take next step
            pred_posA = A + a*b; % prediction of the A's future position
            pred_posB = B - a*b; % prediction of the A's future position
            pred_diff = pred_posA - pred_posB; % difference vector between predicted future position of A and B
            pred_dist = norm(pred_diff);       
            stepsize = stepsize + normrnd([0,0],var_step*[1,1]);
            stepsize(1) = max(min_stepsize,stepsize(1));
            stepsize(2) = max(min_stepsize,stepsize(2));
            theta = thetavar*[randn(1) + thetamod*sign(imag(pred_diff))*exppdf(abs(imag(pred_diff)), 1), randn(1) + thetamod*sign(imag(pred_diff))*exppdf(abs(imag(pred_diff)), 1)];
            A = A + stepsize(1)*exp(1i*theta(1));
            B = B - stepsize(2)*exp(1i*theta(2));
            
            % plotting %
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A)])
                title("detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if norm(A-B) < 2*r
                hold on
                title("collision")
                hold off
                danger_index = norm(A-B) - 2*r;
                danger_vec(length(danger_vec) + 1) = danger_index;
                encounter_classifier = -2;
                break
            end
        %%%%%%%%%%%%%%%% no detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vectors, and updating status
            danger_index = D - 2*r;
            danger_vec(length(danger_vec) + 1) = danger_index;

            % taking next step
            theta = thetavar*randn([2,1]); %betarnd(alpha0*[1 1],beta0*[1 1]);
            %theta = theta*theta_range - 0.5*theta_range; % step-angle
            stepsize = stepsize + normrnd([0,0],var_step*[1,1]); % .*gamrnd((0.008)^-1*[1 1],0.008*[1 1]); % step-size
            stepsize(1) = max(min_stepsize,stepsize(1));
            stepsize(2) = max(min_stepsize,stepsize(2));
            A = A + stepsize(1)*exp(1i*theta(1));
            B = B - stepsize(2)*exp(1i*theta(2));
            D = norm(A-B);

            % plotting %
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A)])
                title("no detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if D < 2*r                                           % collision has occured before evasive action is taken
                hold on; title("collision"); hold off
                danger_index = D - 2*r;
                
                %%% compute ttc (which will now be negative) 
                danger_FEA(i) = danger_index;
                a1 = A0; a0 = A1;
                b1 = B0; b0 = B1;
                Adiff = a1-a0;
                Bdiff = b1-b0;
                eta = real(conj(Adiff - Bdiff)*(a0-b0))/norm(Adiff-Bdiff)^2;
                mu = (norm(a0 - b0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
                ttc = -eta + [-1 1]*sqrt(eta^2-mu);
                ttc = -max(ttc);
                DAFEA(i,:) = ttc*ones(1,NTTC);
                    
                danger_enc_i(length(danger_enc_i) + 1) = danger_index;
                danger_max(i) = danger_index;
                encounter_classifier = 2;                                  % indicates a collision with no evasive action attempted
            end
    end

    enc_type(i) = encounter_classifier;
    danger_max(i) = min(danger_enc_i);
    
    if encounter_classifier == 1                   % i.e. no evasive action and no collision has occured
        danger_max_nodetec(i) = min(danger_enc_i);
        danger_FEA(i) = Inf;    % Inf indicates that there was no evasive action
    end
    
sum(enc_type==2)   % print to keep track of number of pure collisions
end




%%
DAFEAcol1 = DAFEA(:,1)';
save('danger_CR2million.mat','danger_CR')
save('danger_FEA200.mat','danger_FEA')
save('DAFEA400.mat','DAFEAcol1')
%% collected data
%% parameters
pause_length = 1/4/2/2/2/2/2/2;
plotting = 0;
% parameters for distribution of step-angle
m = 5; % increase m --> less variance in stepsize
alpha0 = m;
beta0 = m;

% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.30*2;
b = 0.2/2;

a_init = 0.1*15;
b_init = 15^-1;

var_step = 0.0001;

% parameters relating to detection probability, pd = exp(-s*(D + d)^p)

d = .60;
s = 0.4;   % decrease to make it harder to detect
p = 1.3;   % increase to make it unlikelier to detect for large d, but more likely for small d
plot(linspace(0,10,100),exp(-s*(linspace(0,10,100)+d).^p))

% Collision Avoidance Parameters: Parameter-Modification Factor = exp( -s*((x + d)/k)^-p )
d_ca1 = 0.01;    % increase to reduce PMF
k_ca1 = 1.4;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca1 = 1.6;    % increases effect of avoidance modifier when 
s_ca1 = 0.2

d_ca2 = 0.1;    % increase to reduce PMF
k_ca2 = 2.6;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca2 = 1.8;    % increases effect of avoidance modifier when 
s_ca2 = 2


thetamod = 15;      % determines how quickly the person can change direction
thetavar = 0.1;
r = 0.18;      % collision radius of each person
xinit = 10;   % indicates how far apart they are at the start
sigma = 0.3;  % variance of initial starting position
N=10000;         % number of encounters
NTTC = 1;  % Number of TTC to sample at first evasive action for each encounter
danger_FEA = zeros(1,N);         % danger index at time of first evasive action
DAFEA = zeros(N, NTTC);          % danger at first moment of evasive action; each row contains NTTC sampled ttc values for one encounter
danger_CR = zeros(1,N);          % danger index at time of conflict resolution
danger_nodetec = zeros(1,N);     % highest index of danger for non-detection encounters
enc_type = zeros(1,N);           % vector containing encounter type of each encounter
compute_ttc = 1;                 % tells the algorithm whether or not you want it to estimate ttc distribution for each encounter
min_stepsize = 0.02;
%%
%%
enc_type_test = enc_type;
danger_FEA_test = danger_FEA;
DAFEA_test = DAFEA;
p_pure = 0.0054 % estimated probability of pure collision
pu = 0.4026;
%% based of 10000 iter
p_pure = 0.0016;
sdev = sqrt( p_pure*(1-p_pure)/10000 )
ci = p_pure + [-1,1]*1.96*sdev
p20 = 0.1990
p18_to_25 = [0.1834    0.1914    0.1990    0.2066    0.2132    0.2200    0.2268    0.2327]
%% empirical estimate based on 20 000 iterations
p_pure = 5.0000e-04
p_inter = 0.0019
p_tot = 0.0024
sdev = sqrt( p_tot*(1-p_tot)/20000 )
ci = p_tot + [-1,1]*1.96*sdev
%%
p_pureest = 2.0256e-04
p_interest = 0.00153947

p_tot_est = 0.0017
