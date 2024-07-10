function H = Hmtrx(r)
% H = HMTRX(r) computes the 6x6 system transformation matrix
% 
% S = Smtrx(r);
% H = [eye(3)     S'
%      zeros(3,3) eye(3) ];       Property: inv(H(r)) = H(-r)
%
% If r = r_g is the vector from CO to CG, the model matrices in CO and CG
% are related by: M_CO = H(r_g)' * M_CG * H(r_g). Generalized position and
% force satisfy: eta_CO = H(r_g)' * eta_CG and tau_CO = H(r_g)' * tau_CG 
%
% Author:    Thor I. Fossen
% Date:      2001-06-14
% Revisions: 

S = Smtrx(r);
H = [eye(3)     S'
     zeros(3,3) eye(3) ];
