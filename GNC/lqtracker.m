function [G1,G2,G3] = lqtracker(A,B,C,Q,R,E)
% LQTRACKER  computes the LQ tracker gain matrices for LTI systems:
%     dx/dt = Ax + Bu + Ew where E is an optionally input for diturbance feedforward
%     [G1,G2]    = lqtracker(A,B,C,Q,R)   returns u = G1*x + G2*yd 
%     [G1,G2,G3] = lqtracker(A,B,C,Q,R,E) returns u = G1*x + G2*yd + G3*w
%
% Author: Thor I. Fossen
% Date: 16th June 2001
% Revisions: 

[K,P,EE] = lqr(A,B,C'*Q*C,R);
G1 = -inv(R)*B'*P;
Temp = inv((A+B*G1)');
G2 = -inv(R)*B'*Temp*C'*Q;

if nargin==6,
   G3 = inv(R)*B'*Temp*P*E;
else
   G3 = NaN;
end
