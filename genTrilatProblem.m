function [x,r,B] = genTrilatProblem(p,m,n,s)

%genTrilatProblem / generates a problem (for testing).
%
%[X,R,B] = genTrilatProblem(P,M,N,S), where X is the 'real' coordinates of the
%point, R are perturbed radii, B are beacon positions. P is the number of
%points. N is the number of beacons. M is the number of dimensions. S is the
%standard deviation when perturbing R.
%
%Linus Narva

	x = rand(m,p)*10;
	B = rand(m,n)*10;
	r = zeros(n,p);

	for i = 1:p
		for j = 1:n
			r(j,i) = sqrt(sum((x(:,i)-B(:,j)).^2));
		end
	end

	%Perturb r.
	r = r + s*randn(n,p);
end
