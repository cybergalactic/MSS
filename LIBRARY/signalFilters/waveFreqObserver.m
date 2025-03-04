function [x_next, omega_e] = waveFreqObserver(x, y, w_f, K_f, h)
% waveFreqObserver estimates the wave encounter frequency omega_e using the
% nonlinear observer of Belleter, Galeazzi and Fossen (2015). This is a
% signal-based approach where the wave-induced disturbances are assumed to have
% a dominating (peak) frequency, which is observed for fully developed seas. 
% The nonlinear observer is given by the differential equations:
%
%   x1_dot = x2
%   x2_dot = w_f^2 * (y - x1) - 2 * w_f * x2
%   theta_w_dot = K(A) * x1 * (x2_dot - theta_w * x1)
%
% The encounter frequency estimate omega_e is implicitely given by 
% theta_w = -omega_e^2.
%
% Inputs:
%   x     - Current state vector (3x1 vector), inital value [0; 0; 0.5^2]
%   y     - Motion measurement with wave-induced disturbance (scalar)
%   w_f   - Filter natural frequency, larger than the encounter frequency w_e 
%   K_f   - Observer gain
%   h     - Sampling time
%
% Outputs:
%   x_next  - Propagated state vector at time k+1 (5x1 vector)
%   omega_e - Estimate of the wave encounter frequency at time k
%
% Example usage: 
%   exWaveFreqObserver.m
%
% Reference:
%   D. J. Belleter, R. Galeazzi and T. I. Fossen (2015). Experimenal Verification
%   of a Globally Exponentially Stable Nonlinear Wave Encounter Frequency
%   Estimator. Ocean Engineering 97(15), 48–56.
%
% Author: Thor I. Fossen
% Date: 2024-12-01
% Revision:

% Constants
omega_min = 0.4; % Minimum encounter frequency that can be estimated w_e  

% Estimate of the wave encounter frequency
omega_e = max( omega_min, sqrt(abs(x(3))) );

% RK4 discretization of observer
x_next = rk4(@observer, h, x(1:3), y, K_f, w_f); 

end

%% Observer differential equations
function x_dot = observer(x, y, K_f, w_f)
% Nonlinear wave frequency observer (Belleter, Galeazzi and Fossen, 2015) where
% the encounter frequency estimate omega_e is given by theta_w = -omega_e^2.
% Differential equations:
%
%   x1_dot = x2
%   x2_dot = w_f^2 * (y - x1) - 2 * w_f * x2
%   theta_w_dot = K_f * x1 * (x2_dot - theta_w * x1)
x_dot = [
    x(2)
    w_f^2 * (y - x(1)) - 2 * w_f * x(2)
    K_f * x(1) * ( w_f^2 * (y - x(1)) - 2 * w_f * x(2) - x(3) * x(1) ) ];

end

