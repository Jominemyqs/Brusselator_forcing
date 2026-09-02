function F = brusselator_frozen_rhs(state, cfg, b)
%BRUSSELATOR_FROZEN_RHS Evaluate the semidiscrete Brusselator vector field.
%   F = BRUSSELATOR_FROZEN_RHS(STATE,CFG,B) returns an N-by-2 array whose
%   columns are du/dt and dv/dt. STATE is the repository-standard N-by-3
%   array [x,u,v]. The finite-difference stencil and Neumann boundary rows
%   exactly match SOLVE_BRUSSELATOR_1D_FORCED.

validateattributes(state, {'numeric'}, {'2d','real','finite'});
if size(state,2) ~= 3 || size(state,1) < 3
    error('brusselator_frozen_rhs:InvalidState', ...
        'state must be an N-by-3 [x,u,v] array with N >= 3.');
end
validateattributes(b, {'numeric'}, {'scalar','real','finite'});

x = state(:,1);
if any(diff(x) <= 0)
    error('brusselator_frozen_rhs:InvalidGrid', ...
        'The spatial grid must be strictly increasing.');
end
dx = diff(x);
if max(abs(dx-dx(1))) > 1e-10 * max(1,abs(dx(1)))
    error('brusselator_frozen_rhs:NonuniformGrid', ...
        'The method-of-lines solver assumes a uniform grid.');
end

N = size(state,1);
e = ones(N,1);
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) ...
    - sparse(1:N,1:N,e,N,N);
D2 = D2 + D2';
D2(1,2) = 2;
D2(N,N-1) = 2;
D2 = D2 / dx(1)^2;

u = state(:,2);
v = state(:,3);
a = cfg.model.a;
d2 = cfg.model.d2;
d1 = cfg.model.sigma * d2;
F = [d1*D2*u + a - (b+1)*u + u.^2.*v, ...
     d2*D2*v + b*u - u.^2.*v];
end
