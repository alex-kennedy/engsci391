function [varstatus,basicvars,cB,Binv,xB] = fullupdate(~,c,s,r,BinvAs,phase1,varstatus,basicvars,Binv,xB,n)
% Updates the basis representation.
% Input:
%     m = number of constraints
%     c = nx1 cost vector
%     s = index of entering variable
%     r = position in the basis of the leaving variable
%     BinvAs = mx1 Binv*As vector
%     phase1 = boolean, phase1 = true if Phase 1, or false otherwise
%     varstatus = 1xn vector, varstatus(i) = position in basis of variable i,
%     or 0 if variable i is nonbasic
%     basicvars = 1xm vector of indices of basic variables
%     cB = mx1 basic cost vector
%     Binv = mxm basis inverse matrix
%     xB = mx1 basic variable vector
% Output:
%     varstatus = 1xn updated varstatus vector
%     basicvars = 1xm updated basicvars vector
%     cB = mx1 updated basic cost vector
%     Binv = mxm updated basis inverse matrix
%     xB = mx1 updated basic variable vector

% Gauss-Jordan Pivot
xB(r) = xB(r) ./ BinvAs(r);
Binv(r,:) = Binv(r,:) ./ BinvAs(r);

BinvAs(r) = 0;

xB = xB - xB(r) * BinvAs;
Binv = Binv - Binv(r,:) .* BinvAs;

% Book keeping
% Update varstatus and basicvars
if basicvars(r) <= length(varstatus)
    varstatus(basicvars(r)) = 0;
end
if s <= length(varstatus)
    varstatus(s) = 1;
end

basicvars(r) = s;

if phase1
    cB = (basicvars > n) * 1;
else
    cB = c(basicvars(basicvars <= n));
end

end
