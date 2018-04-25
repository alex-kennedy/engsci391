function  [result,z,x,pi] = fullrsm(m,n,c,A,b)
% Solves a linear program using Gauss-Jordon updates
% Assumes standard computational form
% Performs a Phase I procedure starting from an artificial basis
% Input:
%   m,n     = number of constraints and variables
%   c       = nx1 cost vector
%   A       = mxn constraint matrix
%   b       = mx1 rhs vector
% Output:
%   result  = 1 if problem optimal, 0 if infeasible, -1 if unbounded
%   z       = objective function value
%   x       = nx1 solution vector
%   pi      = mx1 dual vector

% For representation error
tol = 1e-8;

% Initialisations
basicvars = n + 1:n + m;
Binv = eye(m);
varstatus = zeros(1,n);

% Move into Phase I
phase1 = true;

[~,minrc,varstatus,basicvars,Binv,xB,~] = partialrsm(m,n,c,A,b,Binv,varstatus,basicvars,phase1);

% Check Feasibility
if any(basicvars > n)
    % Do any basic variables remain in the basis?
    
    if all((xB(basicvars > n) < tol) & (xB(basicvars > n) > -tol))
    
        % If artificial variables are 0, continue to Phase II
        feasible = true;
    
    else
        
        % Could not drive artificial variables from the basis, infeasible
        feasible = false;
    
    end
elseif minrc > tol
    
    % If all artificial variables were driven from the basis, but the
    % minimum reduced cost is still > 0, declare infeasible
    feasible = false;
    
else
    feasible = true;
end

if ~feasible
    % Problem is infeasible. Could not drive all variables from the basis
    % and reduce cost, or, artifical variables remained with non-zero x
    % values.
    result = 0;
    z = [];
    x = [];
    pi = [];
    return
end

% Begin Phase II
phase1 = false;

[result,~,~,basicvars,~,xB,pi] = partialrsm(m,n,c,A,b,Binv,varstatus,basicvars,phase1);

if result == 1
    x = zeros(n,1);
    x(basicvars(basicvars <= n)) = xB(basicvars <= n);

    z = c.' * x;
else
    % Empty values for unbounded
    z = [];
    x = [];
    pi = [];
end

end
