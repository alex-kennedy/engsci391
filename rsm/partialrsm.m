function [result,minrc,varstatus,basicvars,Binv,xB,pi] = partialrsm(m,n,A,b,c,Binv,varstatus,basicvars,phase1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
while true
    pi = getpi(m,n,c,Binv,basicvars,phase1);
    
    [s,minrc] = fullfindEV(m, n, c, A, phase1, varstatus, pi);
    
    BinvAs = Binv * A(:,s);
    xB = Binv * b;
   
    if minrc >= 0
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

    [varstatus,basicvars,cB,Binv,xB] = fullupdate(m,c,s,r,BinvAs,phase1,varstatus,basicvars,Binv,xB,n);
    
%     disp(cB(basicvars))
%     disp(xB)
%     disp(cB * xB)
    
end
end

