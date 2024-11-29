function [xf_next, y] = notchFilter(xf, u, w_0, zeta, h)
% notchFilter is compatible with MATLAB and GNU Octave (www.octave.org).
% This function implements a 2nd-order notch filter or a cascaded double
% 2nd-order notch filter using the RK4 method. The filter transfer function is:
%
% h_notch(s) = (s^2 + 2*zeta*w_0*s+ w_0^2) / (s^2 + 2*w_0*s + w_0^2)
%
% Inputs:
% xf          - Current state (2x1 for single notch, 4x1 for double notch)
% u           - Input signal to be filtered 
% w_0         - Notch frequency
% zeta        - Damping factor
% h           - Sampling time
%
% Outputs:
% xf_next     - Propagated filter state at time k+1
% y           - Filtered output at time k+1
%
% Author: Thor I. Fossen
% Date: 2024-11-27
% Revision: 

xf(:);

% Define 2nd-order filter state-space matrices
A2 = [0 1; -w_0^2 -2*zeta*w_0];
B2 = [0; 1];
C2 = [w_0^2 2*zeta*w_0];
D2 = 1;

% Define 4th-order (cascaded) filter state-space matrices
A4 = [A2, zeros(2); zeros(2), A2];
B4 = [B2; B2];
C4 = [C2, zeros(1, 2)] * [A2, zeros(2); zeros(2), eye(2)] + [zeros(1, 2), C2];
D4 = D2;

% Determine system order (2nd-order or 4th-order)
if length(xf) == 2
    % 2nd-order filter
    A = A2;
    B = B2;
    C = C2;
    D = D2;
elseif length(xf) == 4
    % 4th-order filter
    A = A4;
    B = B4;
    C = C4;
    D = D4;
else
    error('State vector length must be 2 or 4 for 2nd-order or 4th-order filter.');
end

% RK4 implementation for state propagation
k1 = A * xf + B * u;
k2 = A * (xf + 0.5 * h * k1) + B * u;
k3 = A * (xf + 0.5 * h * k2) + B * u;
k4 = A * (xf + h * k3) + B * u;

xf_next = xf + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

% Compute the filtered output
y = C * xf_next + D * u;

end

