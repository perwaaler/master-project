function out = ttc_simulator(A0,B0,r, thetavar, a, b)
stepcounter = 0;
if real(A0)>=real(B0)
    out = "bad initial choice"
end
while real(A0)<real(B0)
            % generate new step
            theta = thetavar*randn([2,1]);
            stepsize = gamrnd(a*[1 1],b*[1 1]); 
            A1 = A0 + stepsize(1)*exp(1i*theta(1));
            B1 = B0 - stepsize(2)*exp(1i*theta(2));
            % solve system to find ttc, if on collision course
            Adiff = A1-A0;
            Bdiff = B1-B0;
            eta = real((Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
            mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
            ttc = -eta + [-1 1]*sqrt(eta^2-mu^2);
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