pause_length = 2^-8;
plotting = 1; % set to one if plots of encounters are wanted

% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.30*2;
b = 0.2/2;

% parameters for the gamma distribution for the initial steps
a_init = 0.08*10;
b_init = 10^-1;

sigma_step = 0.0001;

% parameters relating to detection probability, p(d) = exp(-s*(D + d)^p)
A = 1;     % amplitude; reduce to make pd smaller everywhere
d = .60;   % acts as a focal point; pd for D larger than this are suppressed by larger p, while D smaller than this are increased by larger p
s = 0.4;   % decrease to make it harder to detect; lower to make difference smaller for large and small seperation distances D
p = 1.3;   % This parameter allows you to regulate difference in pd for large and small D; if p is small then pd becomes more uniform
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


thetamod = 1.2;      % determines how quickly the person can change direction
thetavar = 0.1;
r = 0.06;                        % collision radius of each person
xinit = 10;                      % determines how far apart they are at the start
sigma = 0.3;                     % variance of initial starting position
N = 400;                         % number of encounters
NTTC = 3000;                     % Number of TTC to sample at first evasive action for each encounter
danger_FEA = zeros(1,N);         % danger index at time of first evasive action
DAFEA = zeros(N, NTTC);          % Saves danger at first moment of evasive action; each row contains NTTC sampled ttc values for one encounter
danger_CR = zeros(1,N);          % danger index at time of conflict resolution
danger_nodetec = zeros(1,N);     % highest index of danger for non-detection encounters
enc_type = zeros(1,N);           % vector containing encounter type of each encounter
compute_ttc = 0;                 % tells the algorithm whether or not you want it to estimate ttc distribution for each encounter
min_stepsize = 0.02;             % set to make it impossible for a stepsize to be smaller than this