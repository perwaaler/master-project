function out = mod_fac(x, min_val, dillution, exponent, scale)
% out = mod_fac(x, min_val, dillution, exponent, scale)
% out = exp(-scale*((x + min_val)/dillution)^-exponent)
%
% min_val increases the fcn value at x = 0 
%
% dillution to reduces fcn value for all x
%
% increasing exponent makes fcn valur greater for x > dillution - minval
% and smaller for x < dillution-minval
%
% increase scale to increase fcn value for all x
out = exp(-scale*((x + min_val)./dillution).^-exponent);
end