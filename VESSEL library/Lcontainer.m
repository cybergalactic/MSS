function [xdot,U] = Lcontainer(x,ui,U0)
% [xdot,U] = Lcontainer(x,ui,U0) returns the speed U in m/s (optionally) and the 
% time derivative of the state vector: x = [ u v r x y psi p phi delta ]'  using the
% the LINEARIZED model corresponding to the nonlinear model container.m. 
%
% u     = surge velocity          (m/s)
% v     = sway velocity           (m/s)
% r     = yaw velocity            (rad/s)
% x     = position in x-direction (m)
% y     = position in y-direction (m)
% psi   = yaw angle               (rad)
% p     = roll velocity           (rad/s)
% phi   = roll angle              (rad)
% delta = actual rudder angle     (rad)
%
% The inputs are :
%
% Uo     = service speed (optinally. Default speed U0 = 7 m/s
% ui     = commanded rudder angle   (rad)
%
% Reference:  Son og Nomoto (1982). On the Coupled Motion of Steering and 
%             Rolling of a High Speed Container Ship, Naval Architect of Ocean Engineering,
%             20: 73-83. From J.S.N.A. , Japan, Vol. 150, 1981.
% 
% Author:    Thor I. Fossen
% Date:      23rd July 2001
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>


% Check of input and state dimensions
if (length(x) ~= 9),error('x-vector must have dimension 9 !');end

% Check of service speed
if nargin==2, U0=7.0; end
if U0 <=0,error('The ship must have speed greater than zero');end

% Normalization variables
rho = 1025;                 % water density (kg/m^3)
L = 175;                    % length of ship (m)
U = sqrt(U0^2 + x(2)^2);    % ship speed (m/s)

% rudder limitations
delta_max  = 10;             % max rudder angle (deg)
Ddelta_max = 5;              % max rudder rate (deg/s)

% States and inputs
delta_c = ui(1); 

v = x(2);     y = x(5);
r = x(3);   psi = x(6);
p = x(7);   phi = x(8);
nu    = [v r p]';   
eta   = [y psi phi]';
delta = x(9);    
 
% Linear model using nondimensional matrices and states with dimension (see Fossen 2002): 
% TM'inv(T) dv/dt + (U/L) TN'inv(T) v + (U/L)^2 TG'inv(T) eta = (U^2/L) T b' delta

T    = diag([ 1 1/L 1/L]);
Tinv = diag([ 1 L L ]);

M = [ 0.01497    0.0003525     -0.0002205       
      0.0003525  0.000875       0        
     -0.0002205  0              0.0000210 ];

N = [ 0.012035   0.00522    0   
      0.0038436  0.00243   -0.000213   
     -0.000314   0.0000692  0.0000075  ];

G = [ 0 0 0.0000704  
      0 0 0.0001468
      0 0 0.0004966 ];

b = [-0.002578 0.00126  0.0000855 ]';

% Rudder saturation and dynamics
if abs(delta_c) >= delta_max*pi/180,
   delta_c = sign(delta_c)*delta_max*pi/180;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max*pi/180,
   delta_dot = sign(delta_dot)*Ddelta_max*pi/180;
end

% TM'inv(T) dv/dt + (U/L) TN'inv(T) v + (U/L)^2 TG'inv(T) eta = (U^2/L) T b' delta
nudot = inv(T*M*Tinv)*((U^2/L)*T*b*delta-(U/L)*T*N*Tinv*nu-(U/L)^2*T*G*Tinv*eta );

% Dimensional state derivatives xdot = [ u v r x y psi p phi delta ]'
xdot =[  0
         nudot(1)
         nudot(2)
         cos(psi)*U0-sin(psi)*cos(phi)*v
         sin(psi)*U0+cos(psi)*cos(phi)*v
         cos(phi)*r               
         nudot(3)
         p
         delta_dot                ];

