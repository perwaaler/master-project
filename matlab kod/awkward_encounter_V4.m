clear all
pause_length = 2^-10;


% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.30*2;
b = 0.2/2;

% parameters for the gamma distribution for the initial steps
step_mean = 0.13;
step_var = 0.005;
a_init = step_mean^2/step_var; %0.09*10;
b_init = step_var/step_mean;   %10^-1;

% variance of the noise added to each step
sigma_step = 0.0001;

% parameters relating to detection probability, p(d) = exp(-s*(D + d)^p)
Amp = 0.93;     % amplitude; reduce to make pd smaller everywhere
d = 0;          % acts as a focal point; pd for D larger than this are suppressed by larger p, while D smaller than this are increased by larger p
s = 0.5;   % decrease to make it harder to detect; lower to make difference smaller for large and small seperation distances D
p = 2;   % This parameter allows you to regulate difference in pd for large and small D; if p is small then pd becomes more uniform
plot(linspace(0,10,100),Amp*exp(-(s*linspace(0,10,100)).^p))

% Collision Avoidance Parameters: Parameter-Modification Factor = exp( -s*((x + d)/k)^-p )
d_ca1 = 0.01;    % increase to reduce PMF
k_ca1 = 1.4;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca1 = 1.6;    % increases effect of avoidance modifier when 
s_ca1 = 0.2;


thetamod_k = 0.3;               % determines how quickly the person can change direction after detecting the other person
thetamod_s = 1.6;                % makes modification of behaviour more dramatic for very small ydiff, but less dramatic for large ydiff
thetamod_p = 3;                  % amplifies difference in reaction between small and large ydiff
thetavar = 0.023;                % determins the amount of random variability of the step-angle at each iteration.
thetavar0 = 0.015;               % variance of initial theta
r = 0.06;                         % collision radius of each person
xinit = 6;                      % determines how far apart they are at the start
sigma = 0.3;                     % variance of initial starting position
N =500;                         % number of encounters
NTTC = 200;                      % Number of TTC to sample at first evasive action for each encounter
danger_FEA = zeros(1,N);         % danger index at time of first evasive action
DAFEA = ones(N, NTTC)*-10;          % Saves danger at first moment of evasive action; each row contains NTTC sampled ttc values
danger_CR = zeros(1,N);          % danger index at time of conflict resolution
danger_nodetec = zeros(1,N);     % highest index of danger for non-detection encounters (i.e. type 1 encounters)
enc_type = zeros(1,N);           % vector containing encounter type of each encounter
plotting = 0;                    % set to one if plots of encounters are wanted
compute_ttc =1;                 % tells the algorithm whether or not you want to estimate ttc distribution for each encounter
min_stepsize = 0.02;             % set to make it impossible for a stepsize to be smaller than this
stability_fac = 0.04;
dangermax_EA = nan*ones(N,100);
for i=1:N
i
% initiation of walks
EA_index = 1;
first_detection = 0;                                                 % They have not yet detected eachother
stepsize = min_stepsize*[1 1] + [gamrnd(a_init*[1 1],b_init)];        % determines average speed
theta = normrnd([0,0],thetavar0); 
A_real = -xinit;
B_real = xinit;
A_im = normrnd(0,sigma);
B_im = normrnd(0,sigma);
A0 = A_real+1i*A_im;
B0 = B_real+1i*B_im;
% data collection variables
detection_status = 0; 
encounter_classifier = 1; % 1 --> (no detection, no crash), -1 --> (detection, no crash) 2 --> (no detection, crash), -2 --> (detection, crash)
danger_vec = [];          % collects the danger index associated with each timestep
detec_vec = [];           % tracks detection status, 1 --> detection at timeframe i, 0 --> no detection at timeframe i

    while real(A0) < xinit && real(B0) > -xinit
        D = norm(A0-B0);
        Dim = norm(imag(A0-B0));
        pd = Amp*exp(-(s*D)^p)*exp(-(2.5*Dim)^3); % detection probability
        %%%%%%%%%%%%%%%% detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  (detection_status==1 | rand(1) < pd) && real(A0) < real(B0)
            % saving information in vectors, and updating status
            if first_detection == 0 & compute_ttc == 1 % a loop that collects information at first moment of evasive action
                for k=1:NTTC
                    ttc = ttc_simulator_double_momentum(A0,B0,stepsize,theta,min_stepsize,sigma_step,r, thetavar,stability_fac);
                    DAFEA(i,k) = ttc;
                end
                danger_FEA(i) = norm(A0-B0) - 2*r;
                dangermax_EA(i,EA_index) = norm(A0-B0) - 2*r;
                EA_index = EA_index + 1;
            end
            
            first_detection = 1;                            % set to 1 so that above loop only runs on first moment of detection
            detection_status = 1;
            encounter_classifier = -1;
            detec_vec(length(detec_vec) + 1) = 1;
            danger_index = D - 2*r;
            danger_vec(length(danger_vec) + 1) = danger_index;

            % take next step
            pred_posA = A0 + a*b; % prediction of the A's future position
            pred_posB = B0 - a*b; % prediction of the B's future position
            pred_diff = pred_posA - pred_posB; % difference vector between predicted future position of A and B
            pred_dist = norm(pred_diff);       % distance between predicted positions at time i+1  
            alpha_modifier = exp(-s_ca1*((pred_dist + d_ca1)/k_ca1)^-p_ca1);
            stepsize = stepsize + normrnd([0,0],sigma_step); 
            stepsize(1) = max(min_stepsize,stepsize(1));
            stepsize(2) = max(min_stepsize,stepsize(2));
            theta = normrnd(theta,thetavar) + normrnd(-stability_fac*theta,thetavar) + exp(-0.2*pred_dist)*thetamod_k*sign(imag(pred_diff))*exppdf(abs(imag(thetamod_s*pred_diff))^thetamod_p, 1);
            A1 = A0 + stepsize(1)*exp(1i*theta(1));
            B1 = B0 - stepsize(2)*exp(1i*theta(2));
            dangermax_EA(i,EA_index) = norm(A1-B1) - 2*r;
            EA_index = EA_index + 1;
            
            % plotting %
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)])
                title("detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B1)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if norm(A1-B1) < 2*r                % this part is entered, then a colliosion with evasive action has occured
                hold on
                title("collision")
                hold off
                detec_vec(length(detec_vec) + 1) = 1;
                danger_index = norm(A1-B1) - 2*r;
                danger_vec(length(danger_vec) + 1) = danger_index;
                encounter_classifier = -2;             % set encounter to type crash with attempted evasive action
                break
            end
%%%%%%%%%%%%%%%%%%%%%%%% no detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vectors, and updating status
            detec_vec(length(detec_vec) + 1) = 0;
            
            danger_index = D - 2*r;
            danger_vec(length(danger_vec) + 1) = danger_index;

            % taking next step
            theta = theta + normrnd(-stability_fac*theta,thetavar);
            stepsize = stepsize + normrnd([0 0],sigma_step);
            stepsize(1) = max(min_stepsize,stepsize(1));
            stepsize(2) = max(min_stepsize,stepsize(2));
            A1 = A0 + stepsize(1)*exp(1i*theta(1));
            B1 = B0 - stepsize(2)*exp(1i*theta(2));
            D = norm(A1-B1);

            % plotting %
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)])
                title("no detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B1)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if norm(A1-B1) < 2*r % collision has occured before detection
                hold on
                title("collision")
                hold off
                detec_vec(length(detec_vec) + 1) = 1;
                danger_index = norm(A1-B1) - 2*r;
                    danger_FEA(i) = norm(A1-B1) - 2*r;
                    a1 = A0; a0 = A1;
                    b1 = B0; b0 = B1;
                    Adiff = a1-a0;
                    Bdiff = b1-b0;
                    eta = real(conj(Adiff - Bdiff)*(a0-b0))/norm(Adiff-Bdiff)^2;
                    mu = (norm(a0 - b0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
                    ttc = -eta + [-1 1]*sqrt(eta^2-mu);
                    ttc = -max(ttc);
                    DAFEA(i,:) = ttc*ones(1,NTTC);
                danger_vec(length(danger_vec) + 1) = danger_index; 
                encounter_classifier = 2;                        % indicates a collision with no evase action attempted
                break
            end
        end
        A0 = A1;
        B0 = B1;
    end

    enc_type(i) = encounter_classifier;
    detec_dangerind = danger_vec(find(detec_vec)); % contains danger index for all times with detection status
    if encounter_classifier == 1     % i.e. no evasive action and no collision
        danger_nodetec(i) = min(danger_vec);
        danger_FEA(i) = 1000;    % 1000 indicates that there was no evasive action
    end
    if encounter_classifier == 2
        %danger_CR(i) = detec_dangerind(end);    
    end
sum(enc_type==2)
end

danger_FEA = danger_FEA(find(danger_FEA<1000)); % remove all elements corresponding to encounters with no EA
DAFEA(find(enc_type==1),:) = [];                % remove all elements corresponding to encounters with no EA
dangermax_EA(find(enc_type==1),:) = [];