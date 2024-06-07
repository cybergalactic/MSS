function [xdot, U, M, D] = supply(x, tau)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [xdot, U, M, D] = supply(x, tau) returns the time derivative 
% xdot = A * x + B * tau of the state vector: x = [ x y psi u v r ]' and 
% the speed U = sqrt( u^2 + v^2 ) for a supply vessel length L = 76 m. 
% The 3x3 mass matrix M and 3x3 damping matrix D are optionally outputs, 
% which can be used for control design.
%
% The model is only valid for staionkeeping (dynamic positioning) and 
% low speed, that is U < 3 m/s. The states are:
%
%   u:   surge velocity                    (m/s)     
%   v:   sway velocity                     (m/s)
%   r:   yaw velocity                      (rad/s)
%   x    position in the x-direction       (m)
%   y:   position in the y-direction       (m)
%   psi: yaw angle                         (rad)
%
% The control forces and moment are given by tau = [X, Y, N]'.
%
% Example usage:
%   [~,~,M,D] = supply()      : Return the 3x3 mass ans damping matrices 
%   [xdot, U] = supply(x,tau) : Return xdot and U
%   xdot = supply(x,tau)      : Return xdot and U
%
% Reference: 
%   T. I. Fossen, S. I. Sagatun and A. J. Sørensen (1996). Identification 
%   of Dynamically Positioned Ships. Journal of Control Engineering 
%   Practice CEP-4(3):369-376
%
% Author:     Thor I. Fossen
% Date:       12 July 2002
% Revisions:  
%   24 Feb 2004 - Included missing mass in the Bis transformation
%   12 Oct 2011 - Corrected T and Tinv, which were switched 
%   27 May 2019 - Added U as ouput
%   31 May 2019 - Included the rotation matrix in yaw
%   22 Mar 2023 - Corrected wrong assignmnet of states

if nargin == 0, x = zeros(6,1); tau = zeros(3,1);  end 

% Normalization variables
L    =  76.2;           % Length of ship (m)
g    =  9.81;           % Acceleration of gravity (m/s^2)
mass = 6000e3;          % Mass (kg)

T    = diag([1 1 1/L]);
Tinv = diag([1 1 L]);

% Model matricses
Mbis = [1.1274         0          0
             0    1.8902    -0.0744
             0   -0.0744     0.1278];

Dbis = [0.0358        0        0
             0        0.1183  -0.0124
             0       -0.0041   0.0308];
 
 
 % Check of input and state dimensions
 if (length(x)  ~= 6),error('x-vector must have dimension 6!');end
 if (length(tau) ~= 3),error('u-vector must have dimension 3!');end
 
 psi = x(3);
 R = [ cos(psi) -sin(psi) 0
       sin(psi)  cos(psi) 0
             0         0  1 ];
  
 M = mass * Tinv^2 * (T * Mbis * Tinv);
 D = mass * Tinv^2 * (sqrt(g/L) * T * Dbis * Tinv);
 
 A = [ zeros(3,3)         R
       zeros(3,3) -inv(M)*D ];
 
 B = [ zeros(3,3)
       inv(M) ];
 
 % Dimensional state derivative
 xdot = A * x + B * tau;
 
 % Speed
  U = sqrt( x(4)^2 + x(5)^2 );
 
