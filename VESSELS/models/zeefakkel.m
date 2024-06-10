function [psi_dot, r_dot,delta_dot] = zeefakkel(r,U,delta,delta_c,d_r)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [rdot] = zeefakkel(r,U,delta,delta_c,d_r) returns the yaw acceleration, 
% yaw rate, and rudder angle of the Norrbin (1963) nonlinear autopilot model
%
%                  psi_dot = r                     Yaw kinematics
%  T r_dot + n3 r^3 + n1 r = K delta + d_r         Norrbin (1963) model
%                dot_delta = f(delta, delta_c)     Rudder dynamics
% 
% for the recreational craft Zeefakkel (length 45 m and beam 8 m). The craft 
% is controlled by controllable pitch propeller and a rudder with rudder 
% dynamics. The model parameters K, T and, n3 are speed dependent, while n = 1 
% (course stable ship). The parameters are interpolated for varying speeds U.
%
% Inputs: 
%   r:       yaw rate                  (rad/s)     
%   U:       speed                     (m/s)
%   delta:   rudder angle              (rad)
%   delta_c: commanded rudder angle    (rad)
%   d_r:     yaw moment disturbance    (optional)    
%
% Reference: 
%   J. Van Amerongen (1982). Adaptive Steering of Ships – A Model
%   Reference Approach to Improved Maneuvering and Economical Course 
%   Keeping. PhD thesis. Delft University of Technology, Netherlands.
%
% Author:     Thor I. Fossen
% Date:       2020-06-19
% Revisions:  

if (nargin == 4), d_r = 0; end
if (U < 1 || U > 7), error('U should be between 1-7 m/s'); end

% Model parameters
n1 = 1; 
n3 = 0.4;
delta_max  = deg2rad(30);       % Max rudder angle (rad)
Ddelta_max = deg2rad(10);       % Max rudder rate (rad/s)

% Zeefakkel (Van Amerongen 1982) 
data = [...                              data = [ U K T ]
    2       0.15    33   
    2.6     0.19    33
    3.6     0.29    33
    4       0.37    33
    5       0.50    31
    6.2     0.83    43 ];

% Interpolate to find K and T as a function of U
K = interp1(data(:,1), data(:,2), U, 'linear', 'extrap');
T = interp1(data(:,1), data(:,3), U, 'linear', 'extrap');

% Rudder saturation and dynamics
delta_c = sat(delta_c, delta_max);

delta_dot = delta_c - delta;
delta_dot = sat(delta_dot, Ddelta_max);

% Yaw dynamics
psi_dot = r;
r_dot = (1/T) * (K * delta + d_r - n3 * r^3 - n1 * r ); 

