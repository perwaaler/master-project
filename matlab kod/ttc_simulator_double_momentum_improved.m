function out = ttc_simulator(A0,B0,stepsize0,theta0,min_stepsize,sigma_step, r, thetavar,stability_fac)
% out = ttc_simulator(A0,B0,stepsize0,min_stepsize,sigma_step, r, thetavar)
% This function generates a random walk starting at position A0 and
% B0 respectively, and then computes Time To Collision. TTC is set to 1000
% if collision never occurs.
% A0 and B0 given current positions, and theta describe current direction.

stepcounter = 0;
ttc = 1i;
D = norm(A0-B0);

% if they have already collided, run:
if D < 2*r
    % get pre-collision positions based on arguments
    B1 = B0 + stepsize0(2)*exp(1i*theta0(2));
    A1 = A0 - stepsize0(1)*exp(1i*theta0(1));
    Adiff = A1-A0;
    Bdiff = B1-B0;
    eta = real(conj(Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
    mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
    ttc = -eta + [-1 1]*sqrt(eta^2-mu);
    ttc = -max(ttc); % ttc will now be negative, since time must be reversed to find time of collision
    out = ttc;
else
    % generate walk
    while real(A0) < real(B0)
                D = norm(A0-B0);
                % generate new step
                theta0 = theta0 + normrnd(-stability_fac*theta0,thetavar);
                stepsize0 = stepsize0 + normrnd([0,0],sigma_step);
                stepsize0(1) = max(min_stepsize, stepsize0(1));
                stepsize0(2) = max(min_stepsize, stepsize0(2));
                A1 = A0 + stepsize0(1)*exp(1i*theta0(1));
                B1 = B0 - stepsize0(2)*exp(1i*theta0(2));         

                % collision has occured
                if norm(A1-B1)<2*r
                    Adiff = A1-A0;
                    Bdiff = B1-B0;
                    eta = real(conj(Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
                    mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
                    ttc = -eta + [-1 1]*sqrt(eta^2-mu);
                    ttc = min(ttc);
                    out = stepcounter + ttc;
                    break
                else
                    stepcounter = stepcounter + 1;
                    A0 = A1;
                    B0 = B1;
                end
    end
    
    %if no collision has occured, set:
    if ttc==1i
        out = Inf;
    end
end
end