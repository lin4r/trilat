function [x,code,n,xTrace,alphaTrace]=gn(x0,resfun,params, ...
		convTol,maxIter,minAlpha,armijo,backrate)

%GN \ Globally convergent Gauss-Newton implementation based on
%'Armijo backtracking' line search.
%
%
%
%X=GN(X0,RESFUN,PARAMS), returns the optimum value X. X0 is the
%initial guess. RESFUN is the problem-independent residual (and
%Jacobian of residual) function. PARAMS is a cell array of
%additional parameters to RESFUN (ther than the iterates).
%
%X=GN(X0,RESFUN,PARAMS,CONVTOL,MAXITER,MINALPHA,ARMIJO,BACKRATE),
%also specifies constants to be used (if uspecified then default
%values are used, see DEFAULT VALUES). CONVTOL is the tolerance
%for convergence (see CONVERGENCE TEST). MAXITER is the maximum
%number of iterations allowed. MINALPHA is the minimum backtrack
%step allowed. ARMIJO is the Armijo constant. BACKRATE is the
%update in the backtracked step size when the Armijo condition
%fails.
%
%[X,CODE]=... also returns an error code CODE. 0 on
%convergence, -1 on failure to converge (in MAXITER), -2 if the
%Armijo condition failed.
%
%[X,CODE,N]=... also returns the number of iterations executed
%N.
%
%[X,CODE,N,XTRACE]... also returns the iteration trace XTRACE.
%The X corresponding to each iteration is stored in the columns
%of XTRACE (ordered first to last).
%
%[X,CODE,N,XTRACE,ALPHATRACE]... also returns the trace of line
%search step sizes in a vector ALPHATRACE (ordered first to
%last).
%
%
%
%RESIDUAL FUNCTION INTERFACE
%The residual function RESFUN must have the following body:
%[R,J] = RESFUN(X, ...), where R is the residual and J the
%Jacobian of the residual. The iterate X must be the first
%parameter. Any number of additional parameters is allowed.
%
%
%
%DEFAULT VALUES FOR CONSTANTS
%If the constants are omitted, then the following default values
%will be used.
%
%	convTol = 1e-8
%	maxIter = 50
%	armijo = 1e-4
%	minAlpha = 1e-3
%	backrate = 0.5
%
%
%
%CONVERGENCE TEST
%Convergence is based on the angle between the residual and the
%tangent plane. The cosine of the angle theta is given by
%
%	cos(theta) = || Js || / || r ||
%
%When the cosine is smaller than convTol convergence has
%occurred.
%
%Furthermore the 1 times convTol is added to the denominator
%so that 'perfect' data can be handled (otherwise a zero
%residual would cause division by zero).
%
%The convergence criteria used is then:
%
%	||J*p||<=convTol*(1+||r||)
%
%Linus Narva

	%Use default values if the parameters was not specified.
	if nargin==3
		convTol = 1e-8;
		maxIter = 50;
		armijo = 1e-4;
		minAlpha = 1e-3;
		backrate = 0.5;
	end

	x = x0; %x0 is the first iterate.

	%Allocate the traces.
	xTrace = zeros(length(x),0);
	alphaTrace = [];

	%The convergence test.
	convTest = @(J,r,p,convTol) norm(J*p) ...
		<=convTol*(1+norm(r));

	%The Gauss-Newton step rule. Solve J'*J*p=J*(-r) for p.
	gnStep = @(J,r) J'*J\-J'*r;

	% Get the residual and the Jacobian.
	[r,J] = resfun(x,params{:});

	p = gnStep(J,r);

	alpha = 1; %Initialize so that the loop test succeeds.
	n = 0; %Start at 0 iterations.

	%Until Convergence OR line search condition fails, or
	%maxIter is superseded.
	while ~convTest(J,r,p,convTol) ...
			&& (alpha~=0) ...
			&& (n<maxIter)

		%Find line search step length.
		[alpha,flg] = backtrack(x,resfun,params,p ...
			,armijo,minAlpha,backrate);

		if nargout>4
			alphaTrace = [alphaTrace;alpha];
		end

		%Update x.
		x = x+alpha*p;

		%If requested store the x iteration in xTrace.
		if nargout>3
			xTrace = [xTrace,x];
		end

		%Count iteration.
		n = n+1;

		%Calculate convergence data for next iteration.

		%Get the residual and the Jacobian.
		[r,J] = resfun(x,params{:});

		p = gnStep(J,r);
	end

	%Check for errors.
	if alpha==0
		code = flg;
	elseif (n>=maxIter) && ~convTest(J,r,p,convTol)
		code = -1;
	else
		code = 0;
	end
end

function [alpha,flg] = backtrack(x,resfun,params,p,armijo ...
		,alphaMin,backrate)

%BACKTRACK / Backtracks until the Armijo condition is satisfied.
%
%ALPHA=backtrack(X,RESFUN,PARAMS,P,ARMIJO,ALPHAMIN,BACKRATE),
%returns the step size determined by Armijo backtracking, or 0
%on failure. X is the current iterative value. RESFUN is the
%residual function (see GN). PARAMS is a cell array of
%additional parameters to resfun (other than X). P is the search
%direction. ARMIJO is the Armijo constant. ALPHAMIN is the
%minimum allowed step size. BACKRATE is the rate (a number in
%(0,1)) the factor that updates the step size.
%
%[ALPHA,FLG]=... also returns the FLG from ARMIJOCOND.
%
%Linus Narva (c10lna) 2014-05-05

	%Get the residual and Jacobian.
	[r,J] = resfun(x,params{:});

	%Initialize the step size to 1.
	alpha = 1;

	while 1

		% Test the line search condition.
		[sat,flg]=armijoCond(alpha,x,p,resfun,params ...
			,armijo);
		if sat
			break;
		end;

		% Update step size.
		alpha = backrate*alpha;

		% If alpha is too small, then fail.
		if alpha<alphaMin
			alpha = 0;
			break;
		end
	end
end

function [sat,flg] = armijoCond(alpha,x,p,resfun,params,armijo)

%CONDARMIJO / Tests if an update satisfies the Armijo
%condition.
%
%SAT=CONDARMIJO(ALPHA,X,P,RESFUN,PARAMS,ARMIJO) returns true
%(non-zero) if and only if the Armijo condition was satisfied.
%ALPHA is the step size, X the current iterative value. P the
%(line) search direction. RESFUN the residual function (see GN).
%PARAMS are parameters to RESFUN other than x. ARMIJO is the
%Armijo constant.
%
%[SAT,FLG]=... also returns an error flag (-2 if the condition
%failed), 0 on success.
%
%Linus Narva (c10lna) 2014-05-05

	[r,J] = resfun(x,params{:});

	% Calculate the linear model function of the update.
	F = 0.5*r'*r; % F(x)
	gradF = J'*r; % gradient(F(x))
	m = F+armijo*alpha*gradF'*p; % Linear model for F(t);

	% Calculate the real value of the objective function
	%after the update.
	t = x+alpha*p;
	r_ = resfun(t,params{:});
	F_ = 0.5*r_'*r_; % F(t)

	% Test the condition.
	sat = F_<m;

	% Set the value of flg.
	if sat
		flg = 0;
	else
		flg = -2;
	end
end
