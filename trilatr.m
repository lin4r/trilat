function [f,J] = trilatr(x,r,B)

%trilatr \ Residual function (and Jacobian) for trilateration
%in Rn space.
%
%[F,J] = trilatr(X,R,B), computes the residual F and Jacobian J
%for the coordinates x and length to beacons r (column vectors).
%B is the matrix that contains the coordinates of the beacons
%(row-wise as column vectors).
%
%PROBLEM DESCRIPTION
%Trilateration is the problem of determining the position x of
%an object, knowing the distances r from the object to some
%reference points, at the positions B.
% This corresponds to finding the cross point of circles
%centred at B and with respective radii r.
% In 2D, three beacons (which do not lie on the same line) are
%needed for the problem to have an unique solution. If only two
%reference points are used the problem will have two solutions.
% Similarly, in 3D, four beacons (which do not lie on the
%same plane) are needed for a unique solution. If only three
%reference points are used, then the problem will have 2
%solutions.
%
%Linus Narva

	[d,~] = size(x);
	[~,n] = size(B);

	%Coordinate-wise differences in position between the
	%object and the beacons
	diff = repmat(x,1,n) - B;

	%The estimated radii aka. distances to beacons.
	rr = sqrt(sum(diff.^2)');

	%The residual is the estimated distance minus the
	%measured distance.
	f = rr - r;

	if nargout > 1
		%The Jacobian.
		J = diff'./repmat(rr,1,d);
	end
end
