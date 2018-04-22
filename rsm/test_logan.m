%% Golbon's test cases from Canvas
% Results verified by two scripts.
clear
fprintf('Running tests...\n\n')

%% infeasible
task = 'Detect an infeasible problem';
test = 1;
c = [   1;  -1; 0;  0;  0   ];
A = [   1   1   1   0   0   ;
        1   -1  0   -1  0   ;
        1   0   0   0   1   ];
b = [ 200; 100; 50 ];
[m, n] = size(A);

result_ans = 0;

fprintf('[Test %d] %s...\n', test, task);
[result,~,~,~] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
fprintf('[Test %d] Passed\n\n', test);

%% optimal solution
task = 'Find an optimal solution';
test = 2;
c = [   1;  1;  0;  0   ];
A = [   1   1   1   0   ;
        1   -1  0   -1  ];
b = [ 200; 100 ];
[m, n] = size(A);

result_ans = 1;
z_ans = 100;
x_ans = [ 100; 0; 100; 0 ];
pi_ans = [ 0; 1 ];

fprintf('[Test %d] %s...\n', test, task);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
assert(isequal(z,z_ans), '[Test %d] Solver returned wrong objective z', test);
assert(isequal(x,x_ans), '[Test %d] Solver returned wrong solution vector x', test);
assert(isequal(pi,pi_ans), '[Test %d] Solver returned wrong dual vector pi', test);
fprintf('[Test %d] Passed\n\n', test);

%% unbounded
task = 'Detect an unbounded problem';
test = 3;
c = [   1;  1;  -1; 0;  0   ];
A = [   1   1   -1  1   0   ;
        1   -1  1   0   -1  ];
b = [   200;    100 ];
[m, n] = size(A);

result_ans = -1;

fprintf('[Test %d] %s...\n', test, task);
[result,~,~,~] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
fprintf('[Test %d] Passed\n\n', test);

% Emma's test cases from Facebook
%% optimal solution
task = 'Find an optimal solution';
test = 4;
c = [  -100;    -150;   0;  0   ];
A = [   1,      2,      1,  0   ;
        3,      1.5,    0,  1   ];
b = [ 40; 48 ];
[m, n] = size(A);

result_ans = 1;
z_ans = -3200;
x_ans = [ 8; 16; 0; 0 ];
pi_ans = [ -66.6667; -11.1111 ];    % to 4dp

fprintf('[Test %d] %s...\n', test, task);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
assert(z-z_ans<1e-2, '[Test %d] Solver returned wrong objective z', test);
assert(max(abs(x-x_ans))<1e-2, '[Test %d] Solver returned wrong solution vector x', test);
assert(max(abs(pi-pi_ans))<1e-2, '[Test %d] Solver returned wrong dual vector pi', test);
fprintf('[Test %d] Passed\n\n', test);

%% unbounded
task = 'Detect an unbounded problem';
test = 5;
c = [   0;  -2; -1; 0;  0;  0   ];
A = [   1,  -1, 0,  1,  0,  0   ;
        -2, 1,  0,  0,  1,  0   ;
        0,  1,  -2, 0,  0,  1   ];
b = [ 5; 3; 5 ];
[m, n] = size(A);

result_ans = -1;

fprintf('[Test %d] %s...\n', test, task);
[result,~,~,~] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
fprintf('[Test %d] Passed\n\n', test);

%% infeasible
task = 'Detect an infeasible problem';
test = 6;
c = [   1;  1;  0;  0;  0   ];
A = [   1,  0,  -1, 0,  0   ;
        0,  1,  0,  -1, 0   ;
        1,  1,  0,  0,  1   ];
b = [ 6; 6; 11 ];
[m, n] = size(A);

result_ans = 0;

fprintf('[Test %d] %s...\n', test, task);
[result,~,~,~] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
fprintf('[Test %d] Passed\n\n', test);

%% artificial var left in basis at end of Phase I
task = 'Find an optimal solution with degenerate starting basis in Phase II';
test = 7;
c = [   -1; 0;  0   ];
A = [   1   2   1   ;
        1   2   1   ];
b = [ 1; 1 ];
[m, n] = size(A);

result_ans = 1;
z_ans = -1;
x_ans = [ 1; 0; 0 ];
% pi_ans = [ -1; 0 ];

fprintf('[Test %d] %s...\n', test, task);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(isequal(result,result_ans), '[Test %d] Solver returned wrong result', test);
assert(isequal(z,z_ans), '[Test %d] Solver returned wrong objective z', test);
assert(isequal(x,x_ans), '[Test %d] Solver returned wrong solution vector x', test);
% assert(max(abs(pi-pi_ans))<1e-2, '[Test %d] Solver returned wrong dual vector pi', test);
fprintf('[Test %d] Passed\n\n', test);

fprintf('All tests passed. Chocolate for you.\n');

%%   RSM examples
% How to use: Run (F5)
%
% Author: Minsang Kim (Ryan & Dani's examples added)

clear;
clc;

tol=1e-6;

%% example 1
A = [1 -1 1 1 0 0; 1 1 0 0 -1 0; 0 0 1  0 0 1];
c = [1;2;3;0;0;0];
b = [4;0;6];
[m,n] = size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(abs(z-0)<tol);

%% example 2
A= [2 1 1 1 0; 
    0 1 -1 0 1];
b =[6 3]';
c = [-3 -2 1 0 0]';
[m,n] = size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(abs(z+10.5)<tol);

%% example 3
 A = [1 2 2 1 0 0;
     2 1 2 0 1 0;
     2 2 1 0 0 1];
 b = [20 20 20]';
 c = [-10 -12 -12 0 0 0]';
[m,n] = size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(abs(z+136)<tol);

%% example 4

% class example
A = [3 2 1 2 1 0 0;
   1 1 1 1 0 1 0;
   4 3 3 4 0 0 1];
c = [-19 -13 -12 -17 0 0 0]';
b = [225 117 420]';
[m,n] = size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(abs(z+1827)<tol);

%% example 5
A = [1 1 0 1 1 0 0 1; 0 1 1 1 1 1 0 0; 1 0 1 1 1 0 1 1; 1 1 0 0 0 1 1 1;0 1 1 1 0 1 0 0; 1 1 0 1 0 1 0 1;1 0 0 0 1 0 1 1;0 1 1 1 0 0 0 0; 0 0 0 0 1 0 1 1; 1 1 1 1 0 -1 0 1; 1 1 0 1 1 1 0 1; 0 1 0 0 0 1 0 1; 1 0 0 0 0 0 -1 1; 0 1 1 1 0 0 1 0; 0 0 1 0 1 1 0 0; 1 0 0 1 1 1 1 0];
c = [-5 -6 -3 -7 3 2 -1 -5]';
b = [4 3 6 7 5 6 4 8 6 9 12 4 3 5 7 11]';
A = [A eye(16)];
c = [c; zeros(16,1)];
[m,n] = size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(abs(z+28)<tol);

% for this example you must have Binary.csv file on your working directory
%example 6
% A = csvread('Binary.csv');
% A = [A -eye(302)];
% b = 2*ones(302,1);
% c = [ones(302,1) ;zeros(302,1)];
% [m,n] = size(A);
% [result,z,x,pi] = fullrsm(m,n,c,A,b);
% assert(abs(z-10.488731487433204)<tol);

%% example 9
 A = [-1	1	0	1	-1	0	0	1
0	3	1	1	1	2	0	0
1	0	-1	-1	1	0	1	1
1	1	0	-1	-1	1	0	2
0	1	1	2	1	1	1	0
1	3	1	1	-1	2	-1	1
1	0	0	1	-1	-1	1	1
0	-1	1	-1	0	1	1	-1
-1	0	1	0	1	0	1	1
1	1	1	4	1	-1	0	1
1	1	0	1	-1	-1	3	1
0	3	1	0	0	1	0	1
1	-1	5	1	1	3	1	1
0	1	0	1	4	0	0	0
-1	0	-1	0	1	-1	0	1
1	0	0	-1	-1	-1	1	1];
b = [3 9 6 10 5 6 8 3 6 9 12 4 7 5 6 5]';
c = [-5	-6	-3	-7	3	2	-1	-5]';
A = [A eye(16)];
c = [c; zeros(16,1)];
[m,n] = size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(abs(z+35)<tol);

%% Cory & Dani's examples
%%
c= [1; 2; 3; 0; 0];
A= [[1, -1, 1, 1, 0]; [1, 1, 0, 0, -1];[0, 0, 1, 0,  0]];
[m n]=size(A);
b=[4;0;6];
[result,z,x,pi] = fullrsm(m,n,c,A,b); % test case 1
assert(abs(z-22)<tol); % tests for optimality
%%
A=[[1,-2,-1,0,0];[0,3,4,1,0];[0,6,1,0,1]];
c=[2;3;2;0;0];
b=[4;5;7];
[m n]=size(A);

[result,z,x,pi] = fullrsm(m,n,c,A,b); % test case 2
assert(abs(z-8)<tol); % tests for optimality
%%

A=[[1,1,2,1,0,0];[2,0,3,0,1,0];[2,1,3,0,0,1]];
b=[4;5;7];
c=[-3;-2;-4;0;0;0];
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b); % test case 3
 assert(abs(z+21/2)<tol); % tests for optimality
%%
A=[[3,2,1,1,0];[2,5,3,0,1]];
b=[10;15];
c=[-2,-3,-4,0,0]';
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
%%
c=[-2,-3,-4]';
A=[[3,2,1];[2,5,3]];
b=[10;15];
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert((z+130/7)<tol); 
 %%
b=[-1,-1]'; %Testing infeasibility
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(result==0);



%%  Testing degeneracy in phase1/phase2 transition
A=[[1,0,0,0,0];[0,1,0,0,0];[0,0,1,1,1]];
c=[2;3;0;0;0];
b=[10;4;0];
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert((z-32)<tol); 

% Additions by Dani
%%
% Testing unboundedness
A = [1 -1 -1 0; 1 1 0 -1];
b = [1;2];
c = [-1;-1;0;0];
m = 2; n = 4;
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(result==-1);
%%
% Node-arc incidence matrix example -will always have artificial variable
% at value 0
A = [1 1 1 1 0 0 0 0 0 0 0 0;...
    0 0 0 0 1 1 1 1 0 0 0 0;...
    0 0 0 0 0 0 0 0 1 1 1 1;...
    1 0 0 0 1 0 0 0 1 0 0 0;...
    0 1 0 0 0 1 0 0 0 1 0 0;...
    0 0 1 0 0 0 1 0 0 0 1 0;...
    0 0 0 1 0 0 0 1 0 0 0 1];
b = [3;7;5;4;4;6;1];
c = [4;5;5;3;3;2;6;1;1;2;5;4];
m = 7; n = 12;
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert((z-45)<tol); 

disp('Congrats, you''ve passed the test.');

%% test.m
tol=1e-6;
%%
c= [1; 2; 3; 0; 0];
A= [[1, -1, 1, 1, 0]; [1, 1, 0, 0, -1];[0, 0, 1, 0,  0]];
[m n]=size(A);
b=[4;0;6];
[result,z,x,pi] = fullrsm(m,n,c,A,b); % test case 1
assert(abs(z-22)<tol); % tests for optimality
%%
A=[[1,-2,-1,0,0];[0,3,4,1,0];[0,6,1,0,1]];
c=[2;3;2;0;0];
b=[4;5;7];
[m n]=size(A);

[result,z,x,pi] = fullrsm(m,n,c,A,b); % test case 2
assert(abs(z-8)<tol); % tests for optimality
%%

A=[[1,1,2,1,0,0];[2,0,3,0,1,0];[2,1,3,0,0,1]];
b=[4;5;7];
c=[-3;-2;-4;0;0;0];
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b); % test case 3
 assert(abs(z+21/2)<tol); % tests for optimality
%%
A=[[3,2,1,1,0];[2,5,3,0,1]];
b=[10;15];
c=[-2,-3,-4,0,0]';
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);

c=[-2,-3,-4]';
A=[[3,2,1];[2,5,3]];
b=[10;15];
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert((z+130/7)<tol); 
%% 
b=[-1,-1]'; %Testing infeasibility
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(result==0);

%%

% Testing degeneracy in phase1/phase2 transition
A=[[1,0,0,0,0];[0,1,0,0,0];[0,0,1,1,1]];
c=[2;3;0;0;0];
b=[10;4;0];
[m n]=size(A);
[result,z,x,pi] = fullrsm(m,n,c,A,b)

%%
% Additions by Dani

% Testing unboundedness
A = [1 -1 -1 0; 1 1 0 -1];
b = [1;2];
c = [-1;-1;0;0];
m = 2; n = 4;
[result,z,x,pi] = fullrsm(m,n,c,A,b);
assert(result==-1);
%%
% Node-arc incidence matrix example -will always have artificial variable
% at value 0
A = [1 1 1 1 0 0 0 0 0 0 0 0;...
    0 0 0 0 1 1 1 1 0 0 0 0;...
    0 0 0 0 0 0 0 0 1 1 1 1;...
    1 0 0 0 1 0 0 0 1 0 0 0;...
    0 1 0 0 0 1 0 0 0 1 0 0;...
    0 0 1 0 0 0 1 0 0 0 1 0;...
    0 0 0 1 0 0 0 1 0 0 0 1];
b = [3;7;5;4;4;6;1];
c = [4;5;5;3;3;2;6;1;1;2;5;4];
m = 7; n = 12;
[result,z,x,pi] = fullrsm(m,n,c,A,b);
%%
disp('congrats, you passed');