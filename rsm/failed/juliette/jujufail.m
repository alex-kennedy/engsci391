%%Juliette Test case failure, note error using * inner matrix dimensions
%%must agree
load tester.mat

[result,z,x,pi] = fullrsm(m,n,c,A,b);
z_true
z
xopt
x

assert(all(equal( z_true, z )));
disp('you get em')