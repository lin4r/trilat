%Test Script
%
%Generates a trilateration problem and displays relative errors for solutions
%using the linearization strategy and using a crude initial guess.
%
%Linus Narva

p = 5; m = 3; n = 4; s = 0;

[X,r,B] = genTrilatProblem(p,m,n,s);

X

x0 = X+0.2*randn(m,p);

xIterative = trilat(r,B,x0)
xLinearized = trilat(r,B)

% Euclidean norm of columns.
normcols = @(A) sqrt(sum(A.^2));

% Relative errors.
err = @(x) normcols(x - X)./normcols(X);

errIterative = err(xIterative)
errLinearized = err(xLinearized)
