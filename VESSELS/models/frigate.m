function [psi_dot, r_dot, delta_dot] = frigate(r, U, delta, delta_c, d_r)
% Compatibel with MATLAB and the free software GNU Octave (www.octave.org).
% [rdot] = frigate(r,U,delta,delta_c,d_r) returns the yaw acceleration, 
% yaw rate and rudder angle of the Norrbin (1963) nonlinear autopilot model
%
%                  psi_dot = r                     - yaw kinematics
%  T r_dot + n3 r^3 + n1 r = K delta + d_r         - Norrbin (1963) model
%                dot_delta = f(delta, delta_c)     - rudder dynamics
% 
% for a frigate (Length 100 m). The ship is controlled by controllable
% pitch propeller and a rudder with rudder dynamics. The model parameters 
% K, Tm and n3 are speed dependent, while n = 1 (course stable ship). 
% The parametersare interpolated for varying speeds U.
%
% Inputs: 
%   r:       Yaw rate (rad/s)     
%   U:       Speed (m/s)
%   delta :  Rudder angle (rad)
%   delta_c: Commanded rudder angle (rad)
%   d_r:     (OPTIONAL) Yaw moment disturbance (Nm)    
%
% Reference: 
%   J. Van Amerongen (1982). Adaptive Steering of Ships – A Model Reference 
%   Approach to Improved Maneuvering and Economical Course Keeping. 
%   PhD thesis. Delft University of Technology, Netherlands.
%
% Author:     Thor I. Fossen
% Date:       2020-06-20
% Revisions:  

if (nargin == 4), d_r = 0; end
if (U < 5 || U > 12), error('U should be between 5-12 m/s'); end

% Model parameters
n1 = 1;
delta_max  = deg2rad(30);   % max rudder angle (rad)
Ddelta_max = deg2rad(10);   % max rudder rate (rad/s)
T_delta = 1.0;              % rudder time constant (s)

% Frigate (Van Amerongen 1982) 
data = [...                        data = [ U K T n3]
    6     0.08    20    0.4   
    9     0.18    27    0.6
    12    0.23    21    0.3 ];

% Interpolate to find K and T as a function of U
K  = interp1(data(:,1),data(:,2),U,'linear','extrap');
T  = interp1(data(:,1),data(:,3),U,'linear','extrap');
n3 = interp1(data(:,1),data(:,4),U,'linear','extrap');

% Rudder saturation and dynamics
delta_c = sat(delta_c, delta_max);

delta_dot = (delta_c - delta) / T_delta;
delta_dot = sat(delta_dot, Ddelta_max);

% Yaw dynamics
psi_dot = r;
r_dot = (1/T) * (K * delta + d_r - n3 * r^3 - n1 * r ); 