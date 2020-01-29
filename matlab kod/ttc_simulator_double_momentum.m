function out = ttc_simulator(A0,B0,stepsize0,theta0,min_stepsize,sigma_step, r, thetavar,stability_fac)
% out = ttc_simulator(A0,B0,stepsize0,min_stepsize,sigma_step, r, thetavar)
% This function generates a randomw walk that starting at position A= and
% B0 respectively, and then computes Time To Collision. TTC is set to 1000
% if collition never occurs.

stepcounter = 0;
if real(A0)>=real(B0)
    out = "bad initial choice"
end
while real(A0)<real(B0)
            % generate new step
            theta0 = theta0 + normrnd(-stability_fac*theta0,thetavar);
            stepsize0 = stepsize0 + normrnd([0,0],sigma_step); % .*gamrnd((0.008)^-1*[1 1],0.008*[1 1]); % step-size
            stepsize0(1) = max(min_stepsize, stepsize0(1));
            stepsize0(2) = max(min_stepsize, stepsize0(2));
            A1 = A0 + stepsize0(1)*exp(1i*theta0(1));
            B1 = B0 - stepsize0(2)*exp(1i*theta0(2));            % solve system to find ttc, if on collision course
            Adiff = A1-A0;
            Bdiff = B1-B0;
            eta = real(conj(Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
            mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
            ttc = -eta + [-1 1]*sqrt(eta^2-mu);
            if imag(ttc)==0
                ttc = min(ttc);
                break
            else
                stepcounter = stepcounter+1;
            end
            A0 = A1;
            B0 = B1;
            % plotting
%                 xlim([-10,10])
%                 ylim([-4,4])
%                 plot([r*cos(linspace(0,2*pi,50))+real(A0)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A0)])
%                 title("no detection")
%                 hold on
%                 plot([r*cos(linspace(0,2*pi,50))+real(B0)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B0)])
%                 xlim([-10,10])
%                 ylim([-4,4])
%                 hold off
%                 pause(1)
end
if imag(ttc)==0
    out = ttc + stepcounter;
else
    out = 1000; % 100 indicates that no collision occured
end
end