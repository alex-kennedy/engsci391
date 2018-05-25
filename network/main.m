% EngSci 391 Assignement 2
% Semester 1, 2018
% Author: Alex Kennedy

offerPrice = [45; 60; 50; 50; 55; 65; 40; 30; 30];
generatorCapacity = [1000; 430; 400; 1085; 387; 640; 1750; 800; 885];
lineCapacity = [500; 600; 500; 800; 600; 1000; 300; 200; 800; 300; 850];
demand = [1952; 722; 60; 284; 855; 0; 1078; 225; 617; 0];

[A,b,c] = getProblem(offerPrice,generatorCapacity,lineCapacity,demand);
[objective,cities,generators,duals] = solve(A,b,c)

function [A,b,c] = getProblem(offerPrice,generatorCapacity,lineCapacity,demand)
% Determines A, b, and c for the simplified electricity network, as defined
% in standard computational form. 

    % Define the NZ electricity network in terms of city connections and 
    % associated generators
    cityArcs = [
        1, 2;   % A - H
        2, 3;   % H - NP
        2, 4;   % H - N
        3, 5;   % NP - W
        4, 5;   % N - W
        5, 6;   % W - B
        6, 7;   % B - C
        6, 8;   % B - D    
        7, 8;   % C - D
        9, 8;   % T - D
        10, 9   % M - T
    ];
    generatorArcs = [1; 1; 1; 2; 3; 4; 6; 8; 10];

    [nCityArcs, ~] = size(cityArcs);

    [nCities, ~] = size(demand);
    [nGenerators, ~] = size(generatorCapacity);
    
    % Initialise empty A matrix
    A = zeros(nCities, 0);

    % Connect cities (in each direction)
    b = demand;
    c = zeros(2 * nCityArcs, 1);
    for i = 1:nCityArcs
        A = connectNodes(A, cityArcs(i, 1), cityArcs(i, 2));
    end 
    for i = 1:nCityArcs
        A = connectNodes(A, cityArcs(i, 2), cityArcs(i, 1));
    end 
    
    % Add generator arcs
    c = [c; offerPrice];
    newArcs = zeros(nCities, nGenerators);
    
    for i = 1:nGenerators
        newArcs(generatorArcs(i), i) = 1;
    end
    
    A = [A, newArcs];
    
    % Add upper limits on generator supply
    [A,b,c] = limitArcs(A,b,c,23:31,generatorCapacity);
    
    % Add upper limits on city arcs
    [A,b,c] = limitArcs(A,b,c,1:11,lineCapacity);
    [A,b,c] = limitArcs(A,b,c,12:22,lineCapacity);
    
    % Apply Kirchoff's Laws
    [~,n] = size(A);
    b = [b; 0; 0];
    
    newRows = zeros(2, n);
    
    newRows(1, [3,5,13,15]) = 1;
    newRows(1, [14,16,2,4]) = -1;
    
    newRows(2, [19,7]) = 1;
    newRows(2, [8,18]) = -1;
    newRows(2, 20) = -0.4;
    newRows(2, 9) = 0.4;
    
    A = [A; newRows];
    
    % Power loss between W and B
    [m,n] = size(A);
    A(:,[6, 17]) = 0; % disconnect existing arcs
    arcs = zeros(m, 4);
    
    arcs(5,1) = -1;    % W ->
    arcs(6,1) = 0.95;  % B
    
    arcs(6,2) = -1;    % B ->
    arcs(5,2) = 0.95;  % W
    
    arcs(5,3) = -1;    % W ->
    arcs(6,3) = 0.85;  % B
    
    arcs(6,4) = -1;    % B ->
    arcs(5,4) = 0.85;  % W
    
    A = [A, arcs];
    c = [c; zeros(4, 1)];
    
    [A,b,c] = limitArcs(A,b,c,n+1:n+4,[500;500;500;500]);
    
end

function [A] = connectNodes(A,from,to)
% Adds arcs to the A matrix to connect the from node to the to node
% Args:
%   A:      mxn existing A matrix
%   from:   node index of arc tail (<=m)
%   to:     node index of arc head (<=m)
% Returns:
%   A:      adjusted A matrix

    [m, ~] = size(A);

    arc = zeros(m, 1);
    
    arc(from) = -1;
    arc(to) = 1;

    A = [A, arc];
    
end


function [A,b,c] = limitArcs(A,b,c,arcsToLimit,limits)
% Constrains arc values to some limits
% Args: 
%   A:             mxn A matrix
%   b:             mx1 right hand side vector
%   c:             nx1 cost vector
%   arcsToLimit:   list of arcs to limit
%   limits:        list of arc limit values (same size as arcsToLimit
% Returns:
%   A:             adjusted A matrix
%   b:             adjusted right hand side vector
%   c:             adjusted costs vector

    [m,n] = size(A);
    nLimits = length(limits);
    
    A = [A,zeros(m,nLimits)];
    
    bottomSection = zeros(nLimits,n);
    bottomSection(:,arcsToLimit) = eye(nLimits);
    bottomSection = [bottomSection, eye(nLimits)];
    A = [A; bottomSection];
    
    b = [b; limits];
    c = [c; zeros(nLimits,1)];
    
end


function [objective,cities,generators,duals] = solve(A,b,c)
% Solves the electricity network problem and returns the solutions
% Args:
%   A:  mxn A matrix
%   b   mx1 right hand side vector
%   c   nx1 cost vector
% Returns:
%   objective:      objective value of solution
%   cities:         values for network flows between cities
%   generators:     generator utilisation
%   duals:          cost increase per unit increase in demand at each city

[m,n] = size(A);
[~,~,~,pi] = fullrsm(m,n,c,A,b);

options = optimoptions('linprog','Algorithm','dual-simplex');
[~,xopt,objective] = evalc('linprog(c,[],[],A,b,zeros(n,1),[],options);');

cities = xopt(1:11) - xopt(12:22);
generators = xopt(23:31);
duals = pi(1:10);

end