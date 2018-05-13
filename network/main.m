offerPrice = [45; 60; 50; 50; 55; 65; 40; 30; 30];
generatorCapacity = [1000; 430; 400; 1085; 387; 640; 1750; 800; 885];
lineCapacity = [500; 600; 500; 800; 600; 1000; 300; 200; 800; 300; 850];
demand = [1952; 722; 60; 284; 855; 0; 1078; 225; 617; 0];

[A,b,c] = getA(offerPrice,generatorCapacity,lineCapacity,demand);

[m,n] = size(A);
[result,z,x,~] = fullrsm(m,n,c,A,b);

function [A,b,c] = getA(offerPrice,generatorCapacity,lineCapacity,demand)
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
        8, 9;   % T - D
        9, 10   % M - T
    ];
    generatorArcs = [1; 1; 1; 2; 3; 4; 6; 8; 10];

    [nCityArcs, ~] = size(cityArcs);

    [nCities, ~] = size(demand);
    [nGenerators, ~] = size(generatorCapacity);

    % Initialise A, b, c with root node
    n = nCities + nGenerators;
    A = zeros(n, 1);
    A(1,1) = 1;

    c = 0;

    % Connect cities
    b = demand;
    extraZeros = zeros(2 * nCityArcs, 1);
    c = [c; extraZeros];
    for i = 1:nCityArcs
        A = connectNodes(A, cityArcs(i, 1), cityArcs(i, 2), true);
    end 

    % Connect generators
    b = [b; -generatorCapacity];
    c = [c; offerPrice];
    for i = 1:nGenerators
        A = connectNodes(A, nCities + i, generatorArcs(i), false);
    end

    % Add dummy demand node and connect everything to it
    dummyDemand = sum(generatorCapacity) - sum(demand);
    [A,b,c] = addDummyDemand(A,b,c,dummyDemand);

end

function [A] = connectNodes(A,from,to,bothWays)
    [m, ~] = size(A);

    arc = zeros(m, 1);
    
    arc(from) = -1;
    arc(to) = 1;

    A = [A, arc];

    if bothWays
        A = connectNodes(A,to,from,false);
    end
end

function [A,b,c] = addDummyDemand(A,b,c,dummyDemand)
    [n, m] = size(A);
    b = [b; dummyDemand];
    c = [c; zeros(n, 1)];

    A = [A, -eye(n)];
    A = [A; zeros(1, m), ones(1, n)];
end
