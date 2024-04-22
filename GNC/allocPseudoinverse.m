function u = allocPseudoinverse(K,T,W,tau)
% allocPseudoinverse is compatible with MATLAB and GNU Octave (www.octave.org). 
% The function u = allocPseudoinverse(K,T,W,tau) performs unconstrained 
% control allocation by distributing a generalized force vector to inputs 
% while minimizing force.
%
% The relationship tau = T * K * u, where 'tau' is a generalized force vector
% (dimension n), is used to distribute these forces to the input vector 'u'
% (dimension r, where r >= n). The function minimizes the force f = K * u.
% The unconstrained solution is (Fossen 2021, Section 11.2.2)
% 
%   u = inv(K) * inv(W) * T' * inv(T * inv(W) * T')
% 
% The optimal solution exists if the matrix product T * T' is non-singular. 
% The control inputs can also be quadratic function u = abs(n) * n.
%
% Inputs:
%   K: rxr diagonal matrix of force coeffisients 
%   T: nxr constant thruster configuration matrix
%   W: rxr diagonal matrix weighting (prizing) the different control forces 
%      f = K * u
%
% Outputs:
%   u: control inputs
%
% The thruster configuration matrix T can be computed using thrConfig.m. 
% For instance, a ship with one tunnel thruster and one main propller with
% thruster locations l_x = [lx_1, lx_2] and l_y = [ly_1, ly_2] gives
%
%    T = thrConfig( {'T', 'M'}, l_x, l_y)
%
% Author:    Thor I. Fossen
% Date:      2001-11-03
% Revisions: 

if det(T * T') == 0

    error('The product T * T''is singular'); 

elseif det(W) == 0 || ~isequal( W, diag(diag(W)) )

    error('W = diag( [w1, w2,...,wn] ) must be a positive diagonal matrix'); 

elseif det(K) == 0 || ~isequal( K, diag(diag(K)) )

    error('K = diag( [k1, k2,...,kr] ) must be a positive diagonal matrix'); 

else

    Winv = diag(1 ./ diag(W));
    Kinv = diag(1 ./ diag(K));
    u = Kinv * Winv * T' * inv( T * Winv * T') * tau;

end

end
