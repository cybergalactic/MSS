function LP = makeLowPass(w_f, h)
% makeLowPass creates a first-order low-pass filter OBJECT implemented using
% a nested function. The internal filter state is stored inside the function
% handle and updated automatically at each call.
%
%     T_f * xdot + x = u    with   w_f = 1 / T_f
%
% Discretized exactly:
%
%     x[k+1] = exp(-w_f h) * x[k] + (1 - exp(-w_f h)) * u[k]
%
% Usage:
%   LP.psi = makeLowPass(w_f, h);
%   psi_lowpss  = LP.psi(eta(6));       % Filtered yaw angle at time k
%
% The internal filter state x is persistent inside the nested function.
%
% Author: Thor I. Fossen
% Date: 202-12-06
% Revisions:


x = zeros(size(w_f(:)));   % Initialize persistent filter state
phi = exp(-h * w_f(:));

LP = @updateFilter;

    function y = updateFilter(u)
        % Ensure u is a column with matching size
        u = u(:); 

         % Output at time k (before propagation)
        y = x;

         % Update internal state for next sample k+1
        x = phi .* x + (1 - phi) .* u(:);
    end

end
