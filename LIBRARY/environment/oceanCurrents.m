function [nu_c, nu_c_dot] = oceanCurrents( ...
    Vc, betaVc, gammaVc, attitude, nu2, model, parameters, t)
% [nu_c, nu_c_dot] = oceanCurrents(Vc,betaVc,gammaVc,attitude,nu2)
% [nu_c, nu_c_dot] = oceanCurrents( ...
%     Vc,betaVc,gammaVc,attitude,nu2,'stochastic',parameters,t)
% computes ocean-current velocity and its time derivative in BODY using
% the angle definitions and current models in Fossen (2027, Section 10.3).
%
% The current is spatially uniform and irrotational. Its NED velocity is
%
%   v_c^n = Vc * [ cos(gammaVc) * cos(betaVc)
%                  cos(gammaVc) * sin(betaVc)
%                  sin(gammaVc) ]
%
% Vc is the total speed. The azimuth betaVc is measured from North towards
% East, in the direction towards which the current flows. The elevation
% gammaVc is positive downwards, as defined in Section 10.3. A horizontal
% NED current (gammaVc = 0) can have a vertical BODY component when the
% craft is rolled or pitched.
%
% With R mapping BODY to NED and nu2 = [p q r]', the outputs are
%
%   nu_c     = [R' * v_c^n; 0; 0; 0]
%   nu_c_dot = [-S(nu2) * nu_c(1:3) + R' * v_c_dot^n; 0; 0; 0]
%
% Inputs:
%   Vc:       Total current speed (m/s), nonnegative.
%   betaVc:   Current azimuth in NED (rad).
%   gammaVc:  Current elevation, positive down (rad), in [-pi/2, pi/2].
%   attitude: Euler angles [phi theta psi] (rad), or a scalar-first unit
%             quaternion [eta eps1 eps2 eps3]. Row or column vector.
%   nu2:      BODY angular velocity [p q r] (rad/s). Row or column vector.
%   model:    Optional 'constant' (default) or 'stochastic'.
%   parameters: Required for 'stochastic', with fields
%       h:     Sampling time (s), positive.
%       mean:  Mean values [Vc; betaVc; gammaVc].
%       T:     Three positive correlation time constants (s).
%       sigma: Three nonnegative diffusion coefficients. Units are
%              (m/s)/sqrt(s) for speed and rad/sqrt(s) for the angles.
%       min:   Lower limits [Vc_min; betaVc_min; gammaVc_min].
%       max:   Upper limits [Vc_max; betaVc_max; gammaVc_max].
%       id:    Optional name for an independent current realization.
%              The default is 'default'.
%     Each vector has three elements. Limits must contain the mean, speed
%     limits must be nonnegative, and elevation limits lie in [-pi/2, pi/2].
%     Azimuth limits use a continuous interval within mean(2) +/- pi;
%     e.g. [170,190] degrees is a valid interval across the +/-180 boundary.
%     Set gammaVc, mean(3), min(3), max(3), and sigma(3) to zero to retain a
%     horizontal stochastic current.
%   t:        Simulation time (s), nonnegative. Required for 'stochastic';
%             pass t(i) from the main simulation loop.
%
% Outputs:
%   nu_c:        6x1 current velocity [u_c v_c w_c 0 0 0]' at this step.
%   nu_c_dot:    6x1 BODY derivative at this step.
%
% Constant current:
%   The input current is constant in NED, so v_c_dot^n = 0. No random
%   numbers are drawn and no persistent state is used.
%
% Stochastic current:
%   For each component X of [Vc; betaVc; gammaVc], Section 10.3.2 uses
%
%       dX = -(X - mean)/T * dt + sigma * dW.
%
%   The inputs Vc, betaVc, gammaVc specify the initial current and remain
%   fixed in successive calls. The function retains the evolving current
%   internally, with a separate state for each parameters.id. Initial
%   values are projected onto the limits; azimuth is represented about its
%   mean using ssa. The first call defines the start of the sampling grid.
%
%   The unbounded Gauss-Markov process is discretized exactly over h and
%   then saturated at min and max. Between sample times, speed, azimuth,
%   and elevation are interpolated linearly, preserving their limits.
%   v_c_dot^n follows by differentiating the spherical velocity expression.
%   This gives a finite-step derivative for the simulated current path;
%   ideal white-noise-driven paths have no classical derivative.
%
%   Repeated calls in the same sample interval reuse the stored samples,
%   including calls at integration stages. The BODY rotation is evaluated
%   anew at each call. Skipped intervals are advanced internally.
%   The state restarts automatically if time moves backwards or the
%   initial values or stochastic parameters change. A new run starting
%   at t = 0 after a completed run therefore needs no clear command.
%   Use rng(seed) before a run when a repeatable realization is required.
%
% Examples (all angles are in radians):
%   attitude = [0; deg2rad(10); deg2rad(20)];
%   nu2 = [0; 0; 0.01];
%
% -------------------------------------------------------------------------
% EXAMPLE 1: 2-D CONSTANT CURRENT (HORIZONTAL IN NED)
% -------------------------------------------------------------------------
%   [nu_c, nu_c_dot] = oceanCurrents( ...
%       0.5, deg2rad(30), 0, attitude, nu2, 'constant');
%   % 'constant' is the default and can be omitted.
%
% -------------------------------------------------------------------------
% EXAMPLE 2: 3-D CONSTANT CURRENT
% -------------------------------------------------------------------------
%   % Total speed 0.5 m/s, azimuth 30 deg, elevation 10 deg downwards:
%   [nu_c, nu_c_dot] = oceanCurrents( ...
%       0.5, deg2rad(30), deg2rad(10), attitude, nu2, 'constant');
%
% -------------------------------------------------------------------------
% EXAMPLE 3: 2-D STOCHASTIC CURRENT (HORIZONTAL IN NED)
% -------------------------------------------------------------------------
%   p2D.h = 0.05;
%   p2D.mean = [0.5; deg2rad(30); 0];
%   p2D.T = [100; 100; 100];
%   p2D.sigma = [0.01; deg2rad(0.1); 0];
%   p2D.min = [0.2; deg2rad(20); 0];
%   p2D.max = [0.8; deg2rad(40); 0];
%   p2D.id = 'horizontal';
%   t = 0:p2D.h:10;
%   for i = 1:nTimeSteps
%       [nu_c, nu_c_dot] = oceanCurrents( ...
%           0.5, deg2rad(30), 0, attitude, nu2, 'stochastic', p2D, t(i));
%   end
%
% -------------------------------------------------------------------------
% EXAMPLE 4: 3-D STOCHASTIC CURRENT
% -------------------------------------------------------------------------
%   p3D.h = 0.05;
%   p3D.mean = [0.5; deg2rad(30); deg2rad(5)];
%   p3D.T = [100; 100; 100];
%   p3D.sigma = [0.01; deg2rad(0.1); deg2rad(0.05)];
%   p3D.min = [0.2; deg2rad(20); deg2rad(-5)];
%   p3D.max = [0.8; deg2rad(40); deg2rad(15)];
%   p3D.id = 'threeD';
%   t = 0:p3D.h:10;
%   for i = 1:nTimeSteps
%       [nu_c, nu_c_dot] = oceanCurrents( ...
%           0.5, deg2rad(30), deg2rad(5), attitude, nu2, ...
%           'stochastic', p3D, t(i));
%   end
%
% -------------------------------------------------------------------------
% EXAMPLE 5: HYBRID INPUT, HORIZONTAL SPEED Vc PLUS VERTICAL COMPONENT vD_c
% -------------------------------------------------------------------------
% Vc is the horizontal speed, betaVc is the azimuth, and vD_c is the 
% positive-down NED component. Inside the craft function use:
%
%   Vc = 0.5;                 % Horizontal speed (m/s).
%   betaVc = deg2rad(30);     % Azimuth (rad).
%   vD_c = 0.1;                % Vertical NED component, positive down (m/s).
%   Vtotal = hypot(Vc, vD_c);
%   gammaVc = atan2(vD_c, Vc);
%   [nu_c, nu_c_dot] = oceanCurrents( ...
%       Vtotal, betaVc, gammaVc, attitude, nu2, 'constant');
%
% This gives v_c^n = [Vc * cos(betaVc); Vc * sin(betaVc); vD_c]
%
% Reference:
%   T. I. Fossen (2027). Handbook of Marine Craft Hydrodynamics and Motion
%       Control, 3rd ed., Wiley. Section 10.3.
%
% Author: Thor I. Fossen
% Date:   2026-11-02

% -------------------------------------------------------------------------
% INPUT VALIDATION
% -------------------------------------------------------------------------
narginchk(5, 8);
if nargin < 6 || isempty(model)
    model = 'constant';
end
if isstring(model) && isscalar(model)
    model = char(model);
end
if ~ischar(model) || size(model,1) ~= 1 || ...
        ~any(strcmp(model, {'constant', 'stochastic'}))
    error('oceanCurrents:InvalidModel', ...
        'model must be ''constant'' or ''stochastic''.');
end

validateattributes(Vc, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'Vc', 1);
validateattributes(betaVc, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'betaVc', 2);
validateattributes(gammaVc, {'numeric'}, ...
    {'real', 'finite', 'scalar', '>=', -pi/2, '<=', pi/2}, ...
    mfilename, 'gammaVc', 3);
validateattributes(attitude, {'numeric'}, ...
    {'real', 'finite', 'vector'}, mfilename, 'attitude', 4);
validateattributes(nu2, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', 3}, mfilename, 'nu2', 5);

% -------------------------------------------------------------------------
% BODY-TO-NED ROTATION MATRIX
% -------------------------------------------------------------------------
attitude = attitude(:);
nu2 = nu2(:);
if numel(attitude) == 3
    R = Rzyx(attitude(1), attitude(2), attitude(3));
elseif numel(attitude) == 4
    R = Rquat(attitude);
else
    error('oceanCurrents:InvalidAttitude', ...
        'attitude must contain three Euler angles or four unit-quaternion components.');
end

% -------------------------------------------------------------------------
% CURRENT MODEL
% -------------------------------------------------------------------------
current = [Vc; betaVc; gammaVc];
currentRate = zeros(3,1);
if strcmp(model, 'stochastic')
    if nargin < 7
        error('oceanCurrents:MissingParameters', ...
            'The stochastic model requires a parameters structure.');
    end
    if nargin < 8
        error('oceanCurrents:MissingTime', ...
            'The stochastic model requires simulation time t.');
    end
    validateattributes(t, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 't', 8);
    parameters = stochasticParameters(parameters);
    [current, currentRate] = stochasticCurrent(current, parameters, t);
end

% -------------------------------------------------------------------------
% BODY-FRAME CURRENT VELOCITY AND DERIVATIVE
% -------------------------------------------------------------------------
[v_c_n, v_c_dot_n] = currentVector(current, currentRate);
v_c_b = R' * v_c_n;
nu_c = [v_c_b; zeros(3,1)];
nu_c_dot = [-Smtrx(nu2) * v_c_b + R' * v_c_dot_n; zeros(3,1)];
end

% -------------------------------------------------------------------------
% LOCAL FUNCTION: NED CURRENT VELOCITY
% -------------------------------------------------------------------------
function [v_c_n, v_c_dot_n] = currentVector(current, currentRate)
Vc = current(1); betaVc = current(2); gammaVc = current(3);
v_c_n = Vc * [cos(gammaVc) * cos(betaVc)
              cos(gammaVc) * sin(betaVc)
              sin(gammaVc)];
% Columns are the derivatives with respect to speed, azimuth, and elevation.
J = [cos(gammaVc)*cos(betaVc), -Vc*cos(gammaVc)*sin(betaVc), -Vc*sin(gammaVc)*cos(betaVc)
     cos(gammaVc)*sin(betaVc),  Vc*cos(gammaVc)*cos(betaVc), -Vc*sin(gammaVc)*sin(betaVc)
     sin(gammaVc),             0,                          Vc*cos(gammaVc)];
v_c_dot_n = J * currentRate;
end

% -------------------------------------------------------------------------
% LOCAL FUNCTION: STOCHASTIC PARAMETER VALIDATION
% -------------------------------------------------------------------------
function parameters = stochasticParameters(input)
required = {'h', 'mean', 'T', 'sigma', 'min', 'max'};
if ~isstruct(input) || ~isscalar(input) || ~all(isfield(input, required))
    error('oceanCurrents:InvalidParameters', ...
        'parameters must contain h, mean, T, sigma, min, and max.');
end
validateattributes(input.h, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'parameters.h');
parameters.h = input.h;
for k = 2:numel(required)
    name = required{k};
    value = input.(name);
    validateattributes(value, {'numeric'}, ...
        {'real', 'finite', 'vector', 'numel', 3}, ...
        mfilename, ['parameters.' name]);
    parameters.(name) = value(:);
end
if any(parameters.T <= 0) || any(parameters.sigma < 0)
    error('oceanCurrents:InvalidParameters', ...
        'Correlation times must be positive and diffusion coefficients nonnegative.');
end
if any(parameters.min > parameters.mean) || ...
        any(parameters.mean > parameters.max) || parameters.min(1) < 0 || ...
        parameters.min(3) < -pi/2 || parameters.max(3) > pi/2 || ...
        parameters.min(2) < parameters.mean(2)-pi || ...
        parameters.max(2) > parameters.mean(2)+pi
    error('oceanCurrents:InvalidLimits', ...
        'Limits must contain the mean and satisfy the speed and angle ranges.');
end
parameters.id = 'default';
if isfield(input, 'id')
    id = input.id;
    if isstring(id) && isscalar(id)
        id = char(id);
    end
    if ~ischar(id) || isempty(id) || size(id,1) ~= 1
        error('oceanCurrents:InvalidId', ...
            'parameters.id must be a nonempty character vector or scalar string.');
    end
    parameters.id = id;
end
end

% -------------------------------------------------------------------------
% LOCAL FUNCTION: PERSISTENT CURRENT REALIZATIONS
% -------------------------------------------------------------------------
function [current, currentRate] = stochasticCurrent(initial, parameters, t)
persistent streams
if isempty(streams)
    streams = struct('id', {}, 'initial', {}, 'parameters', {}, ...
        'origin', {}, 'lastTime', {}, 'sampleIndex', {}, 'current', {}, 'next', {});
end
index = find(strcmp({streams.id}, parameters.id), 1);
restart = isempty(index);
if ~restart
    timeTolerance = 16 * eps(max([1, abs(t), abs(streams(index).lastTime)]));
    restart = t < streams(index).lastTime - timeTolerance || ...
        ~isequal(initial, streams(index).initial) || ...
        ~isequal(parameters, streams(index).parameters);
end

% -------------------------------------------------------------------------
% INITIALIZATION AND AUTOMATIC RESTART
% -------------------------------------------------------------------------
if restart
    if isempty(index)
        index = numel(streams) + 1;
    end
    current = initial;
    current(2) = parameters.mean(2) + ssa(current(2) - parameters.mean(2));
    current = min(parameters.max, max(parameters.min, current));
    stream.id = parameters.id;
    stream.initial = initial;
    stream.parameters = parameters;
    stream.origin = t;
    stream.lastTime = t;
    stream.sampleIndex = 0;
    stream.current = current;
    stream.next = stochasticStep(current, parameters);
else
    stream = streams(index);
end

% -------------------------------------------------------------------------
% ADVANCE SAMPLES AND INTERPOLATE AT THE REQUESTED TIME
% -------------------------------------------------------------------------
samplePosition = (t - stream.origin) / parameters.h;
sampleTolerance = 16 * eps(max(1, abs(samplePosition)));
sampleIndex = floor(samplePosition + sampleTolerance);
while stream.sampleIndex < sampleIndex
    stream.current = stream.next;
    stream.next = stochasticStep(stream.current, parameters);
    stream.sampleIndex = stream.sampleIndex + 1;
end
fraction = min(1, max(0, samplePosition - stream.sampleIndex));
current = (1-fraction) * stream.current + fraction * stream.next;
currentRate = (stream.next - stream.current) / parameters.h;
stream.lastTime = t;
streams(index) = stream;
end

% -------------------------------------------------------------------------
% LOCAL FUNCTION: BOUNDED GAUSS-MARKOV UPDATE
% -------------------------------------------------------------------------
function currentNext = stochasticStep(current, parameters)
decay = exp(-parameters.h ./ parameters.T);
noiseStd = parameters.sigma .* sqrt( ...
    -(parameters.T/2) .* expm1(-2*parameters.h ./ parameters.T));
currentNext = parameters.mean + decay .* (current - parameters.mean) ...
    + noiseStd .* randn(3,1);
currentNext = min(parameters.max, max(parameters.min, currentNext));
end
