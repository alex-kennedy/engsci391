function [r,minratio] = fullfindLV(n,xB,BinvAs,phase1,basicvars)
% Returns the position in the basis of the leaving variable,
% or returns 0 if no leaving variable exists
% Input:
%   m = number of constraints
%   xB = mx1 basic variable vector
%   BinvAs = mx1 vector of Binv*As
%   phase1 = boolean, phase1 = true if Phase 1, or false otherwise
%   basicvars = 1xm vector of indices of basic variables
% Output:
%   r = position in the basis of the leaving variable
%   minratio = minimum ratio from ratio test

% Set a tolerance to avoid rounding errors breaking code
tol = 1e-8;

ratio = xB ./ BinvAs;

% Begin with extended leaving variable criterion
if ~phase1
    positions = transpose(basicvars > n) & BinvAs ~= 0; 
    if any(positions)
        r = find(positions, 1);
        minratio = 0;
        return
    end
end

% Ensure invalid entries never return as the minimum except in the
% unbounded case
ratio(BinvAs <= tol) = Inf;

[minratio,r] = min(ratio);

% If there are multiple minima, ensure the leaving variable is an
% artificial, if present (i.e. the largest index)
if phase1
    r = find(basicvars == max(basicvars(ratio == minratio)));
end

end

