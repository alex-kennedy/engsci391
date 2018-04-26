function [r,minratio] = fullfindLV(n,xB,BinvAs,phase1,basicvars)
% Returns the position in the basis of the leaving variable,
% or returns 0 if no leaving variable exists.
% Author: Alex Kennedy | aken327 | 460783474
% Input:
%   m = number of constraints
%   xB = mx1 basic variable vector
%   BinvAs = mx1 vector of Binv*As
%   phase1 = boolean, phase1 = true if Phase 1, or false otherwise
%   basicvars = 1xm vector of indices of basic variables
% Output:
%   r = position in the basis of the leaving variable
%   minratio = minimum ratio from ratio test

tol = 1e-8;

ratio = xB ./ BinvAs;

% Begin with extended leaving variable criterion
if ~phase1
    positions = transpose(basicvars > n) & (abs(BinvAs) > tol); 
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
% artificial, if present
artificial = basicvars > n;
if any(artificial)
    r_forced = find(ratio <= minratio + tol & artificial.', 1);
    
    % Check if any artificial r could be found
    if r_forced
        r = r_forced;
    end
    
end

end

