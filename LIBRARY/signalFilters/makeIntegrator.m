function I = makeIntegrator(h, x0)
% makeIntegrator creates a discrete-time integrator OBJECT implemented using
% a nested function. The internal state is stored by the returned function
% handle and updated automatically at each call.
%
%     xdot = u
%
% Discretized using forward Euler integration:
%
%     x[k+1] = x[k] + h * u[k]
%
% Usage:
%   I.z = makeIntegrator(h);               % Scalar zero initial state
%   I.z = makeIntegrator(h, zeros(3,1));   % Vector initial state
%   z_int = I.z(error);                    % Integral state at time k
%
% Inputs:
%   h  - Sampling time (s)
%   x0 - Optional scalar or vector initial state; default 0
%
% Output:
%   I  - Stateful integrator function handle
%
% Author: Thor I. Fossen
% Date: 2026-10-03
% Revisions:

if nargin < 2
    x0 = 0;
end

x = x0(:); % Internal integral state

I = @updateIntegrator;

    function y = updateIntegrator(u)
        % Ensure the input is a column with the same size as the state.
        u = u(:);
        if length(u) ~= length(x)
            error('makeIntegrator:InputSize', ...
                'The input and integral state must have the same size.');
        end

        % Return the state at sample k, then propagate it to sample k+1.
        y = x;
        x = x + h * u;
    end

end
