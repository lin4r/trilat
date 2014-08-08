function [x,flg,its] = trilat(r,B,x0)

%trilat / Solves trilateration using Gauss-Newtons method.
%
%X = trilat(R,B). Assuming that there are p points to determine, the space has
%m dimensions and that the number of beacons is n. X is m*p and contains the
%determined coordinates of the p points. R is n*p and contains radii from
%the n beacons to the p points. B is m*n and contains the coordinates of the n
%beacons. To call this version of the function it is required that n > m.
%
%X = trilat(R,B,X0), X0 will be the 'initial guess' of Gauss Newton. This
%version of the method is applicable when n >= m. In the case where n = m,
%trilateration will have two solutions. Then the correct value of the two
%solutions will be obtained if x0 is sufficiently close to it. Of course in the
%general case X0 must be sufficiently close to the solution for if the function
%is to converge to a solution at all.
%
%[X,FLG]=..., also returns the flags from GN, corresponding to
%each point.
%
%[X,FLG,ITS]=..., also returns the number of GN iterations.
%
%
%
%DEPENDENCIES
%
%  trilatlin.m
%  trilatr.m
%  gn.m
%
%Linus Narva

	if nargin <= 2
		x0 = trilatlin(r,B);
	end

	[m,p] = size(x0);

	% Allocate space.
	x = zeros(m,p);
	flg = zeros(p,1);
	its = zeros(p,1);

	for i = 1:p
		ri = r(:,i);
		[xi,flgi,itsi]=gn(x0(:,i),@trilatr,{ri,B});
		x(:,i) = xi;
		flg(i) = flgi;
		its(i) = itsi;
	end
end
