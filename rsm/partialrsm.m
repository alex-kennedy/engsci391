function [result,minrc,varstatus,basicvars,Binv,xB,pi] = partialrsm(m,n,c,A,b,Binv,varstatus,basicvars,phase1)
% Solves an LP given a Binv is already known. 
%   This function factors out common aspects of Phase I and II. It will
%   solve an LP if Binv is already known. For Phase I, this is set to the
%   identity. For Phase II, a Binv is determined from Phase I. 
% Input:
%   m,n         = number of constraints and variables
%   c           = nx1 cost vector
%   A           = mxn constraint matrix
%   b           = mx1 rhs vector
%   Binv        = mxm basis inverse matrix
%   varstatus   = variable position in basis; 0 otherwise
%   basicvars   = 1xm vector of basis indices
%   phase1      = bool, 1 if solving Phase I; 0 otherwise
% Output: 
%   result      = 1 if problem optimal, 0 if infeasible, -1 if unbounded
%   minrc       = minimum reduced cost of last iteration
%   varstatus   = updated varstatus, as above
%   basicvars   = updated basicvars, as above
%   Binv        = updated Binv, as above
%   xB          = nx1 current basis solution
%   pi          = mx1 dual vector

tol = 1e-8;

while true
    pi = getpi(m,n,c,Binv,basicvars,phase1);
    
    [s,minrc] = fullfindEV(m, n, c, A, phase1, varstatus, pi);
    
    BinvAs = Binv * A(:,s);
    xB = Binv * b;
   
    if minrc >= -tol
        % current bfs is optimal, exit
        result = 1;
        break
    end
    
    [r,minratio] = fullfindLV(n,xB,BinvAs,phase1,basicvars);

    if minratio == Inf
        % problem is unbounded
        result = -1;
        break
    end

    [varstatus,basicvars,~,Binv,xB] = fullupdate(m,c,s,r,BinvAs,phase1,varstatus,basicvars,Binv,xB,n);

end
end

