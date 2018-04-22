function [s,minrc] = fullfindEV(~,~,c,A,phase1,varstatus,pi)
% Returns the index of the entering variable and it's reduced cost,
% or returns 0 if no entering variable exists
% Input:
%   n = number of variables
%   c = nx1 cost vector
%   A = mxn constraint matrix
%   varstatus = 1xn vector, varstatus(i) = position in basis of variable i,
%   or 0 if variable i is nonbasic
%   phase1 = boolean, phase1 = true if Phase 1, or false otherwise
% Output:
%   s = index of the entering variable
%   minrc = reduced cost of the entering variable

if phase1
    reducedCosts =  -pi.' * A;
    reducedCosts(:,logical(varstatus)) = Inf;
else
    reducedCosts = c.' - pi.' * A;
    reducedCosts(logical(varstatus)) = Inf;
end

[minrc,s] = min(reducedCosts);

end

