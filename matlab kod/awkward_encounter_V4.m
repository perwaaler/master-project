% simulation of encounters between two vehicles.

clear all
pause_length = 2^-5;


% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.30*2;
b = 0.2/2;

% parameters for the gamma distribution for the initial steps
step_mean = 0.13;
step_var = 0.005;
% parametes as fcn of mean and variance
a_init = step_mean^2/step_var;
b_init = step_var/step_mean;

% variance of the step-size
var_step = 0.0001;

% parameters relating to detection probability dp; dp = A*exp(-s*(D + d)^p)
% where D is distance between vehicles
Amp = 0.93;     % amplitude; reduce to make dp smaller everywhere
d = 0;          % determines where dp is made smaller and larger by p; dp gets larger for D smaller than d, and smaller for D greater than d
s = 0.5;        % increase to make detection less likely
p = 2;          % makes dp larger r smaller depending on whether D>d or D<d. See above comment.
%plot(linspace(0,10,100),Amp*exp(-(s*linspace(0,10,100)).^p))

thetamod_k = 0.3;                       % determines how quickly the person can change direction after detecting the other person
thetamod_s = 1.6;                       % makes modification of behaviour more dramatic for very small ydiff, but less dramatic for large ydiff
thetamod_p = 3;                         % amplifies difference in reaction between small and large ydiff
thetavar = 0.023;                       % determins the amount of random variability of the step-angle at each iteration.
thetavar0 = 0.015;                      % variance of initial theta
r = 0.06;                               % collision radius of each person
xinit = 6;                              % determines how far apart they are at the start
sigma = 0.3;                            % variance of initial starting position
N = 500;                                % number of encounters
NTTC = 200;                             % Number of TTC to sample at first evasive action for each encounter
danger_FEA = nan*zeros(1,N);            % danger index at time of first evasive action
DAFEA = ones(N, NTTC)*-10;              % Saves danger at first moment of evasive action; each row contains NTTC sampled ttc values
danger_max_nodetec = nan*zeros(1,N);    % highest index of danger for non-detection encounters (i.e. type 1 encounters)
danger_max = nan*zeros(1,N)             % collects data from most dangerous timeframe from each encounter
enc_type = zeros(1,N);                  % vector containing encounter type of each encounter
plotting = 0;                           % set to one if plots of encounters are wanted
compute_ttc = 1;                        % tells the algorithm whether or not you want to estimate ttc distribution for each encounter
min_stepsize = 0.02;                    % set to make it impossible for a stepsize to be smaller than this
stability_fac = 0.04;
danger_max_EA = nan*ones(N,1000);       % collects danger measurements at each time frame where there is evasive action

for i=1:N
i
%%% initiation of encounter
EA_index = 1;                                                         % variable used to find most dangerous moment during attempt to avoid collision
first_detection = 0;                                                  % They have not yet detected eachother
stepsize = min_stepsize*[1 1] + [gamrnd(a_init*[1 1],b_init)];        % determines average speed of each vehicle
theta = normrnd([0,0], thetavar0);
%%% initial positions
A_real = -xinit;
B_real = xinit;
A_im = normrnd(0,sigma);
B_im = normrnd(0,sigma);
A0 = A_real + 1i*A_im;
B0 = B_real + 1i*B_im;
%%% data collection variables
detection_status = 0;
encounter_classifier = 1; % tracks status: sign indicates interaction status (+ --> no interaction, - --> interaction), number indicates collision status (1 --> no collision, 2 --> collision)
danger_enc_i = [];        % collects the danger index associated with each timestep

    while real(A0) < xinit && real(B0) > -xinit
        D = norm(A0-B0);
        Dim = norm(imag(A0-B0));                  %
        dp = Amp*exp(-(s*D)^p)*exp(-(2.5*Dim)^3); % detection probability
%%%%%%%%%%%%%%%%%%%%%%%% detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  (detection_status==1 | rand(1) < dp) && real(A0) < real(B0)
            % saving information in vectors, and updating status
            if first_detection == 0 & compute_ttc == 1 % a loop that collects information at first moment of evasive action
                for k=1:NTTC
                    ttc = ttc_simulator_double_momentum(A0,B0,stepsize,theta,min_stepsize,var_step,r, thetavar, stability_fac);
                    DAFEA(i,k) = ttc;
                end
                danger_FEA(i) = norm(A0-B0) - 2*r;
                danger_max_EA(i,EA_index) = norm(A0-B0) - 2*r;
                EA_index = EA_index + 1;
            end

            first_detection = 1;                            % set to 1 so that above loop only runs on first moment of detection
            detection_status = 1;
            encounter_classifier = -1;
            danger_index = D - 2*r;
            danger_enc_i(length(danger_enc_i) + 1) = danger_index;

            % take next step
                pred_posA = A0 + a*b; % prediction of the A's future position
                pred_posB = B0 - a*b; % prediction of the B's future position
                pred_diff = pred_posA - pred_posB; % difference vector between predicted future position of A and B
                pred_dist = norm(pred_diff);       % distance between predicted positions at time i+1
                stepsize = stepsize + normrnd([0,0],var_step);
                stepsize(1) = max(min_stepsize,stepsize(1));
                stepsize(2) = max(min_stepsize,stepsize(2));
                theta = normrnd(theta,thetavar) + normrnd(-stability_fac*theta,thetavar) + exp(-0.2*pred_dist)*thetamod_k*sign(imag(pred_diff))*exppdf(abs(imag(thetamod_s*pred_diff))^thetamod_p, 1);
                A1 = A0 + stepsize(1)*exp(1i*theta(1));
                B1 = B0 - stepsize(2)*exp(1i*theta(2));

            danger_max_EA(i,EA_index) = norm(A1-B1) - 2*r;
            EA_index = EA_index + 1;

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
                danger_index = norm(A1-B1) - 2*r;
                danger_max(i) = danger_index;
                danger_enc_i(length(danger_enc_i) + 1) = danger_index;
                encounter_classifier = -2;             % set encounter to type crash with attempted evasive action
                break
            end
%%%%%%%%%%%%%%%%%%%%%%%% no detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vectors, and updating status

            danger_index = D - 2*r;
            danger_enc_i(length(danger_enc_i) + 1) = danger_index;

            % taking next step
            theta = theta + normrnd(-stability_fac*theta,thetavar);
            stepsize = stepsize + normrnd([0 0],var_step);
            stepsize(1) = max(min_stepsize,stepsize(1));
            stepsize(2) = max(min_stepsize,stepsize(2));
            A1 = A0 + stepsize(1)*exp(1i*theta(1));
            B1 = B0 - stepsize(2)*exp(1i*theta(2));
            D = norm(A1-B1);

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
                break
            end
        end
        A0 = A1;
        B0 = B1;
    end

    enc_type(i) = encounter_classifier;
    danger_max(i) = min(danger_enc_i);

    if encounter_classifier == 1                   % i.e. no evasive action and no collision has occured
        danger_max_nodetec(i) = min(danger_enc_i);
        danger_FEA(i) = Inf;    % Inf indicates that there was no evasive action
    end
sum(enc_type==2)   % print to keep track of number of pure collisions
end

% remove all elements/rows corresponding to encounters with no EA
danger_FEA = danger_FEA(find(danger_FEA<1000));
DAFEA(find(enc_type==1),:) = [];
danger_max_EA(find(enc_type==1),:) = [];

danger_max_EA = min(danger_max_EA'); % find most dangerous moment during attempt to avoid collision
danger_max_nodetec = danger_max_nodetec(find(danger_max_nodetec ~= nan));
