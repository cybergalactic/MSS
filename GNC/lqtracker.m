function [G1,G2] = lqtracker(A,B,C,Q,R)
% LQTRACKER  computes the LQ tracker gain matrices for LTI systems:
%     dx/dt = Ax + Bu
%     [G1,G2] = lqtracker(A,B,C,Q,R)   returns u = G1*x + G2*yd 
%
% Author: Thor I. Fossen
% Date: 16 June 2001
% Revisions: 9 June 2020 T. I. Fossen - removed disturbance FF 

[K,P,EE] = lqr(A,B,C'*Q*C,R);
G1 = -inv(R)*B'*P;
G2 = -inv(R)*B'*inv((A+B*G1)')*C'*Q;

