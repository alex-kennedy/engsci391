function  [result,z,x,pi] = fullrsm(m,n,c,A,b)
% Solves a linear program using Gauss-Jordon updates
% Assumes standard computational form
% Performs a Phase I procedure starting from an artificial basis
% Input:
%   m,n = number of constraints and variables
%   c = nx1 cost vector
%   A = mxn constraint matrix
%   b = mx1 rhs vector
% Output:
%   result = 1 if problem optimal, 0 if infeasible, -1 if unbounded
%   z = objective function value
%   x = nx1 solution vector
%   pi = mx1 dual vector

tol = 1e-8

% Initialisations
basicvars = n + 1:n + m;
Binv = eye(m);
varstatus = zeros(1,n);

% Move into Phase I
phase1 = true;

[~,minrc,varstatus,basicvars,Binv,~,~] = partialrsm(m,n,A,b,c,Binv,varstatus,basicvars,phase1);

if minrc > -tol
    % Problem is infeasible, could not drive all artificial variables from
    % the basis
    result = 0;
    z = 0;
    x = 0;
    pi = 0;
    return
end

% Begin Phase II
phase1 = false;

[result,~,~,basicvars,~,xB,pi] = partialrsm(m,n,A,b,c,Binv,varstatus,basicvars,phase1);

x = zeros(n,1);
x(basicvars(basicvars <= n)) = xB(basicvars <= n);

z = c.' * x;

end
