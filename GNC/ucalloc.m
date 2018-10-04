function u = ucalloc(K,T,W,tau)
% u = ucalloc(K,T,W,tau) unconstrained control allocation. The generalized
%     force vector tau = T*K*u (dim n) is distributed to the input vector u
%     (dim r) where r>=n by minimizing the force f=K*u.
%
%     An unconstrained solution u = inv(K)*inv(W)*T'*inv(T*inv(W)*T')
%     exists if T*T' is non-singular. 
%
%     - K is a diagonal rxr matrix of force coeffisients 
%     - T is a nxr constant configuration matrix. 
%     - W is a rxr positive diagonal matrix weighting (prizing) the 
%       different control forces f = K*u.
%
% Author:    Thor I. Fossen
% Date:      3rd November 2001
% Revisions: 25th September 2002 - function name changed from alloc.m to ucalloc.m

if det(T*T')==0, 
    error('T*T''is singular'); 
elseif det(W)==0,
    error('W must be positive'); 
else
    Winv = inv(W);
    u = inv(K)*Winv*T'*inv(T*Winv*T')*tau;
end


