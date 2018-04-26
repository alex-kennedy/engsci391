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

tol = 1e-8;

s = 0;
minrc = Inf;

if phase1
    
    for i = find(~varstatus)
        rc = -pi.' * A(:,i);
        if rc < minrc + tol
            s = i;
            minrc = rc;
        end
    end
    
else
    
    for i = find(~varstatus)
        rc = c(i) - pi.' * A(:,i);
        if rc < minrc + tol
            s = i;
            minrc = rc;
        end
    end
    
end

end

