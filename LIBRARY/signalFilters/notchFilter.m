function [xf_next, y] = notchFilter(xf, u, w_0, zeta, h, method)
% notchFilter is compatible with MATLAB and GNU Octave (www.octave.org).
% Notch filters can be applied to Inertial Navigation System (INS) measurements
% to remove the dirst-order wave-inuced motions (wave filtering). This function 
% implements a 2nd-order notch filter using either:
%   - A second-order notch filter discretized by uisng the RK4 method.
%   - A direct IIR filter implementation using difference equations.
%
% The second-order notch filter transfer function is:
%
%    h_notch(s) = (s^2 + 2*zeta*w_0*s+ w_0^2) / (s^2 + 2*w_0*s + w_0^2)
%
% Inputs:
%   xf          - Current state (2x1 vector for RK4, 4x1 for IIR)
%   u           - Input signal to be filtered (scalar)
%   w_0         - Notch frequency (from FFT or waveFreqObserver.m)
%   zeta        - Damping factor, typical value for a narrow notch is 0.01 to 0.1
%                 for the second-order filter, while the IIR can use 0.5.
%   h           - Sampling time
%   method      - 'RK4' (default) or 'IIR'
%
% Outputs:
%   xf_next     - Propagated filter state at time k+1
%   y           - Filtered output at time k
%
% Example usage: 
%   exINSwaveFilter.m
%
% Author: Thor I. Fossen
% Date: 2024-11-27
% Revision: 
%   2025-02-11 : Added option for IIR notch filter

if nargin < 6
    method = 'RK4'; % Default method
end

switch method
    case 'RK4'
        xf = xf(:);

        % State-space matrices
        A = [-2*w_0, -w_0^2; 
             1, 0];
        B = [1; 0];
        C = [2*w_0*(zeta-1), 0];
        D = 1;

        % Filtered output at time k 
        y = C * xf + D * u;

        % RK4 implementation for state propagation
        k1 = A * xf + B * u;
        k2 = A * (xf + 0.5 * h * k1) + B * u;
        k3 = A * (xf + 0.5 * h * k2) + B * u;
        k4 = A * (xf + h * k3) + B * u;

        xf_next = xf + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

    case 'IIR'

        % Stability factor
        r = exp(-zeta * w_0 * h);  
        
        % IIR coefficients
        b0 = 1;
        b1 = -2*cos(w_0 * h);
        b2 = 1;
        a0 = 1;
        a1 = -2*r*cos(w_0 * h);
        a2 = r*r;

        % Store coefficients
        b = [b0, b1, b2];
        a = [a0, a1, a2];

        % Extract previous states
        x_prev = xf(1);    % x[n-1]
        x_prev2 = xf(2);   % x[n-2]
        y_prev = xf(3);    % y[n-1]
        y_prev2 = xf(4);   % y[n-2]
        
        % Compute new output 
        y = b(1)*u + b(2)*x_prev + b(3)*x_prev2 - a(2)*y_prev - a(3)*y_prev2;
        
        % Update filter states 
        xf_next = [u; x_prev; y; y_prev];

    otherwise
        error('Invalid method. Choose ''RK4'' or ''IIR''.');
end

end