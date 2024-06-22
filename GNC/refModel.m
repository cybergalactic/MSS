function [x_d,v_d,a_d] = refModel(x_d,v_d,a_d,x_ref,v_max,zeta_d,w_d,h,eulerAngle)
% Poisition, velocity and acceleration reference model based on the method
% by Fossen (2021, Chapter 12.1.1). 
%
% Inputs:  
%   x_d: current desired position at time t_k
%   v_d: current desired velocity at time t_k
%   a_d: current desired velocity at time t_k
%   x_ref: commanded position 
%   v_max: maximum velocity
%   zeta_d: desired relative damping factor
%   w_d: desired natural frequency (rad/s)
%   h: sampling time (s)
%   eulerAngle: 1 if x_d is an Euler angle, 0 else
%
% Outputs:  
%   x_d: desired position at time: t_k+1 = t_k + h
%   v_d: desired velocity at time: t_k+1 = t_k + h
%   a_d: desired velocity at time: t_k+1 = t_k + h
%  
% Author:    Thor I. Fossen
% Date:      2024-02-05

% Use smallest signed angle for Euler angle errors
if eulerAngle == 1
    e_x = ssa( x_d - x_ref );
else
    e_x = x_d - x_ref;
end

% desired jerk (Fossen 2021, Equation 12.10)
a_d_dot = -w_d^3 * e_x - (2*zeta_d + 1) * w_d^2 * v_d...
    - (2*zeta_d + 1) * w_d * a_d;

% Propagation of states: x_d[k+1], v_d[k+1], a_d[k+1]
x_d = x_d + h * v_d;
v_d = v_d + h * a_d;
a_d = a_d + h * a_d_dot;

% Velocity saturation
if abs(v_d) > v_max
    v_d = sign(v_d) * v_max;
end

end


