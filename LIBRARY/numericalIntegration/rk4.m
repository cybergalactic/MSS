function x = rk4(Function, h, x, varargin)
% RK4 integration 
%
% Inputs:
%   Function - handle to the dynamics function: xdot = Function(x, varargin)
%   h - sampling time in seconds
%   x - current state vector x[k]
%   varargin - additional parameters including the current input vector 
%
% Outputs:
%   x - updated state vector x[k+1]
%
% Example usage:
% 
%   Vehicle dynamics: 
%      x = [ u v w p q r x y z phi theta psi ]'
%      xdot = otter(x, n, mp, rp, V_c, beta_c) 
%   RK4:
%      x = rk4(@otter, h, x, n,mp,rp,V_c,beta_c) 
%
% Author:    Thor I. Fossen
% Date:      2024-07-09
% Revisions: 

% Compute k1
k1 = Function(x, varargin{:});

% Compute k2
k2 = Function(x + 0.5 * h * k1, varargin{:});

% Compute k3
k3 = Function(x + 0.5 * h * k2, varargin{:});

% Compute k4
k4 = Function(x + h * k3, varargin{:});

% Compute the weighted average of the slopes
x = x + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

end