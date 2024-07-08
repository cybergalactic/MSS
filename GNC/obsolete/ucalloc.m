function u = ucalloc(K, T, W, tau)
% u = ucalloc(K,T,W,tau) is obsolete and has been replaced by 
% u = allocPseudoinverse(K,T,W,tau) for unconstrained control allocation.
%
% Author:    Thor I. Fossen
% Date:      3 Nov 2001
% Revisions: 
%   25 Mar 2024 - function name changed to allocPseudoinverse.m

disp('u = ucalloc(K, T, W, tau) is obsolete.')
disp('Use u = allocPseudoinverse(K, T, W, tau)')

u = allocPseudoinverse(K, T, W, tau);

end


