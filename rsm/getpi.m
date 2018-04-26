function [pi] = getpi(m,n,c,Binv,basicvars,phase1)
%Gets the vector of duals, pi
%   During Phase I, costs are 1 if corresponding to an artificial variable;
%   0 otherwise. 
%   If any artificial variables remain in the basis, their costs are set to
%   0 to ensure no index out of range errors occur. They are dealt with by
%   the extended leaving variable criterion in Phase II. 
% Author: Alex Kennedy | aken327 | 460783474
% Inputs:
%   m,n = number of constraints and variables
%   c = nx1 cost vector
%   Binv = mxm basis inverse matrix
%   basicvars = 1xm vector of indices of basic variables
%   phase1 = boolean, true if Phase I, or false otherwise
% Outputs: 
%   pi = 1xn array of duals

if phase1
    pi = (((basicvars > n) .* 1) * Binv).';
else
    if any(basicvars > n)
        c = [c; zeros(m, 1)];
    end
    pi = (c(basicvars).' * Binv).';
end

end

