function x = trilatlin(r,B)

%trilatlin / solves trilateration by using a linearization tool.
%
%X = trilatlin(R,B), Cosider the case with n beacons, m dimensions and
%p points to determine. X is m*p and contains the coordinates of each point
%to be solved for. R is n*p and contains distances from beacons to points.
%B is m*n and contains coordinates of beacons.
%
%The position(s) are calculated by solving the linear system A*p = b, for:
%
% A = 
%
%	x2-x1 y2-y1 z2-z1
%	x3-x1 y3-y1 z3-z1
%	...
%	xm-x1 ym-y1 zm-z1
%
%
% p =
%
%	x-x1
%	y-y1
%	z-z1
%
%
% b =
%
%	0.5*(r1^2 - r2^2 + d21^2)
%	0.5*(r1^2 - r3^2 + d31^2)
%	...
%	0.5*(r1^2 - rm^2 + dm1^2)
%
%Where x,y,z are coordinates of the point to determine, xi, yi zi for i from
%2 to m are beacon coordinates. ri's are distances from beacons and dij is
%the distance from beacon i to beacon j.
%
%Linus Narva

	[m,n] = size(B); % Number of dimensions m and number of beacons n.
	[~,q] = size(r); % Number of points to determine q.

	%Choose an arbitratry linearization tool b1.
	b1 = B(:,1);

	%Remove the linearization tool from B.
	Btail = B(:,2:end);

	%Calculate distances from beacons to the linearization beacon.
	D = Btail-repmat(b1,1,n-1);
	%We don't bother taking the squareroot since it's not used.
	dsquared = sum(D.^2);
	dsquared = dsquared'; %As column vector.

	%A is now the transpose of D.
	A = D';

	%Pick the radii corresponding to the linearization tool.
	r1 = r(1,:);

	%Remove the tail of the radii.
	rtail = r(2:end,:);
	b = 0.5*(repmat(r1,n-1,1).^2 - rtail.^2 + repmat(dsquared,1,q));

	%The coordinates can now be obtained from solving the linear system.
	p = A\b;
	x = p+repmat(b1,1,q);

end
