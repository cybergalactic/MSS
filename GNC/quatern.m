function [J,J1,J2] = quatern(q)
% [J,Rq,Tq] = quatern(q) computes the quaternion transformation matrices
%
%  J = [ Rq     0
%           0  Tq ]
%
% where J1 = Rq and J2 = Tq, see Rquat.m and Tquat.m.
%
% Author:    Thor I. Fossen
% Date:      14 Jun 2001
% Revisions: 07 May 2022  - added calls to Rquat and Tquat

J1 = Rquat(q);
J2 = Tquat(q);
 
J = [ J1  zeros(3,3);
      zeros(4,3) J2 ];