function J = jacobian(fun,x,p,h)
%jacobian approximates the jacobian of FUN using central difference.
%
%   J = jacobian(FUN,X,P,H), where FUN is a function name or handle. X is the
%   vector variable do differentiate over, P is a cell array of additional
%   parameters to FUN. H is the scalar size of the differential. If H is a
%   vector, then H(i) is used along X(i).
%
%   J = jacobian(FUN,X), uses P = {} and H = 1e-6.
%
%   J = jacobian(FUn,X,P), uses the same H as above.
%
%Linus Narva

    if nargin < 3
        p = {};
    end

    if nargin < 4
        h = 1e-6
    end

    y0 = feval(fun,x,p{:});

    n = length(x);
    m = length(y0);

    if isscalar(h)
        h = repmat(h,n);
    end

    J = sparse(m,n);

    I = speye(n);

    for i = 1:n
        J(:,i) = (feval(fun,x+h(i)*I(:,i),p{:}) - ...
            feval(fun,x-h(i)*I(:,i),p{:})) / (2*h(i));
    end

end
