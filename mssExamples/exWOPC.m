function exWOPC()
% exWOPC simulates fully actuated and underactuated weather-optimal
% positioning control (WOPC) for marine craft using nonlinear PID control.
% The example is compatible with MATLAB and Octave (www.octave.org).
% The vessel state is x = [x_n y_n psi u v r]' and supply.m is used for
% the vessel dynamics. 
%
% Run exWOPC with no inputs. Select the actuation configuration using
% the flag below: 0 = fully actuated, 1 = underactuated.
% Underactuated operation uses paired main propellers and one bow tunnel
% thruster, which couples sway force and yaw moment.
%
% Actuation is instantaneous and unsaturated. The example computes thrust
% and propeller speeds; it does not estimate shaft power or power savings.
%
% References:
%   Fossen, T. I. (2027). Virtual Anchor Nonlinear PID Control for Weather-
%   Optimal Positioning of Fully Actuated and Underactuated Marine Craft. 
%   Submitted to Ocean Engineering
%
%   Fossen (2027), T. I. "Handbook of Marine Craft Hydrodynamics and Motion
%   Control," 3rd edition, John Wiley & Sons, Ltd., Chichester, UK.
%
%   Fossen, T. I., S. I. Sagatun and A. J. Sørensen (1996). Identification 
%   of Dynamically Positioned Ships. Journal of Control Engineering 
%   Practice CEP-4(3):369-376.
%
% MSS dependencies:
%   supply.m             - 76.2 m supply vessel (Fossen et. al., 1996)
%   Rzyx.m               - Rotation matrix
%   ssa.m                - Smallest signed angle
%
% Author:     Thor I. Fossen
% Date:       2026-08-12
% Revisions:
%   2026-10-18: Revised gains, 180-min horizon, ideal actuation, and
%               physical sway-yaw coupling.

close all;

flag = 0;                        % 0 = fully actuated, 1 = underactuated

% ------------------------------------------------------------------------------
%% SIMULATION PARAMETERS
% ------------------------------------------------------------------------------
T_final = 10800;                 % Final simulation time (s)
h = 0.05;                        % Sampling time (s)
t = (0:h:T_final)';               % Column time vector
nTimeSteps = length(t);

% WOPC GEOMETRY AND ENVIRONMENT
R_d = 50;                          % Desired virtual-circle radius (m)
p_d = [0; 0];                      % Desired North-East position (m)

% Constant environmental load of 100 kN towards west.
w = -1;
w_n = 1e5 * [0; w; 0];           % Constant load in NED coordinates

% WOPC GAINS 
[~,~,M,D] = supply();               % MSS supply-vessel mass matrix

% Nonlinear PID pole-placement initialization using separate surge, sway,
% and yaw bandwidths. The radial task uses omega_u and the angular task uses
% omega_r. The internal sway bandwidth appears only in derivative damping.
% A diagonal derivative gain is used in this case study to avoid
% sway-yaw derivative cross-coupling and permit K_d(2,2) = 0 in the
% underactuated variant:
%   K_p = Omega_task*M_task*Omega_task,
%   K_d = 2*diag(zeta.*omega_n.*diag(M)),
%   Lambda = Omega_task/10, K_i = K_p*Lambda,
%   k_0 = min(diag(Lambda))/2,
%   K_0 = M(1,1)*k_0^2*I_2.
omega_n = [0.05; 0.1; 0.025];     % Bandwidths (rad/s)
zeta = [1.0; 0.5; 1.0];            % [zeta_u,zeta_v,zeta_r]'
T_d = [-1, 0, 0; 0, -1/R_d, -1];
M_task = (T_d * (M \ T_d')) \ eye(2);
Omega_task = diag([omega_n(1),omega_n(3)]);
K_p = Omega_task * M_task * Omega_task;
K_d = 2 * diag(zeta .* omega_n .* diag(M));
if flag == 1
    K_d(2,2) = 0;
end
Lambda = Omega_task / 10;
K_i = K_p * Lambda;              % Explicit equivalent PID gain
k_0 = min(diag(Lambda)) / 2;     % Slower absolute-position anchor
K_0 = M(1,1) * k_0^2 * eye(2);   % Cartesian position stiffness (N/m)

% THRUSTERS AND CONTROL ALLOCATION (FOSSEN 2027, Chapter 12)
%   #1-#2 Bow tunnel thrusters (RPM)
%   #3-#4 Stern tunnel thrusters (RPM)
%   #5-#6 Main propellers (RPM)
% Ideal actuation: instantaneous thrust without saturation.
% Columns 1--4 are transverse tunnel thrusters, while columns 5--6 are
% the main propellers. The coefficients include thrust gains and moment
% arms, and map u = |n|.*n into generalized force tau = B*u.
B = [0,   0,    0,     0,   31.2,   31.2; ...
     3.2, 3.2,  3.2,   3.2,  0,      0; ...
     96,  70.4, -70.4, -96, -249.6, 249.6];

% INITIAL CONDITIONS
eta = [50; -20; 0];               % [x_n,y_n,psi]'
nu = zeros(3,1);                  % [u,v,r]'
p_0 = p_d;                        % Initial virtual-circle center
e_I = zeros(2,1);                 % WOPC integral state
rho_min = 1.0;                    % Guard for the polar singularity (m)

% Hybrid angular-error lift. The boundary tolerance selects the positive
% rotation at the antipodal initial condition; subsequent samples use the
% standard one-dimensional phase-unwrapping recursion.
angularError = 0;
relativeAnglePrevious = 0;
antipodalHysteresis = deg2rad(5);

% ------------------------------------------------------------------------------
%% STORED WOPC RESPONSE
% ------------------------------------------------------------------------------
% Named fields make the saved simulation signals explicit
wopcData.flag = flag;
wopcData.M = M;
wopcData.D = D;
wopcData.loadNED = w_n;
wopcData.headingEquilibrium = atan2(-w_n(2),-w_n(1));
wopcData.Kp = K_p;
wopcData.Ki = K_i;
wopcData.Kd = K_d;
wopcData.Lambda = Lambda;
wopcData.k0 = k_0;
wopcData.K0 = K_0;
wopcData.time = t;
wopcData.nu = zeros(nTimeSteps,3);
wopcData.eta = zeros(nTimeSteps,3);
wopcData.p0 = zeros(nTimeSteps,2);
wopcData.positionError = zeros(nTimeSteps,2);
wopcData.xError = zeros(nTimeSteps,2);
wopcData.integralError = zeros(nTimeSteps,2);
wopcData.sigma = zeros(nTimeSteps,2);
wopcData.n = zeros(nTimeSteps,6);
wopcData.nCommand = zeros(nTimeSteps,6);
wopcData.tauControl = zeros(nTimeSteps,3);
wopcData.tauDisturbance = zeros(nTimeSteps,3);

% ------------------------------------------------------------------------------
%% MAIN LOOP
% ------------------------------------------------------------------------------
for i = 1:nTimeSteps

    % Polar WOPC coordinates. The relative angle is lifted continuously below 
    % to avoid switching at the +/-pi cut
    delta_p = eta(1:2) - p_0;
    rho = max(norm(delta_p),rho_min);
    gamma = atan2(delta_p(2),delta_p(1));
    x_2 = gamma - eta(3);

    % Lift x_2 + pi from the circle to a continuous real-valued error.
    % MSS ssa() implements wrapping to the principal interval. The first
    % sample follows the q_0 initialization in the paper; later samples
    % use the standard one-dimensional phase-unwrapping recursion.
    relativeAngle = x_2 + pi;
    if i == 1
        q0 = ssa(relativeAngle);
        if q0 < -pi + antipodalHysteresis
            angularError = q0 + 2*pi;
        else
            angularError = q0;
        end
    else
        angleIncrement = ssa(relativeAngle - relativeAnglePrevious);
        angularError = angularError + angleIncrement;
    end
    relativeAnglePrevious = relativeAngle;
    x_tilde = [rho - R_d; angularError];

    % U(x) and T(x)
    R_x2 = [cos(x_2), -sin(x_2); sin(x_2), cos(x_2)];
    H = diag([1,rho]);
    U = H \ R_x2';
    T = [U, [0; -1]];

    % Absolute position error is stored for verification and plotting
    positionError = eta(1:2) - p_d;
    R_nb = Rzyx(0,0,eta(3));

    % Composite error and nonlinear PID law with Cartesian position energy:
    %   sigma = x_tilde + Lambda*e_I
    %   g_p = L'*R_z(psi)'*K_0*(p-p_d)
    %   tau_WOPC = -(K_d*nu + T(x)'*K_p*sigma + g_p)
    % Expanding K_p*sigma gives the PID relation K_i = K_p*Lambda
    sigma = x_tilde + Lambda * e_I;
    g_p = [R_nb(1:2,1:2)' * K_0 * positionError; 0];
    generalizedForce = K_d * nu + T' * K_p * sigma + g_p;

    tau_WOPC = -generalizedForce;
    if flag == 1
        % Bow tunnel thruster #1 (x = +30 m) and equally commanded main
        % propellers #5--6. Stern and remaining bow thrusters are inactive.
        % B_eff = [1 0; 0 1; 0 30] maps [total surge thrust; bow force].
        u_command = zeros(6,1);
        u_command(1) = tau_WOPC(3) / B(3,1);
        u_command(5:6) = tau_WOPC(1) / (B(1,5) + B(1,6));
    else
        u_command = B' * ((B * B') \ tau_WOPC);
    end
    n = sign(u_command) .* sqrt(abs(u_command));
    n_c = n;                       % Instantaneous actuator response
    tau = B * u_command;

    % The constant NED load is transformed to BODY coordinates to
    % form the unknown w(psi) in the vessel model
    w_body = R_nb' * w_n;

    % Supply-vessel dynamics, using its constant M and D matrices.
    nu_dot = M \ (tau + w_body - D * nu);
    e_Idot = x_tilde;                                  % Integral error
    % The absolute-position correction removes the translated equilibrium
    % family without introducing a desired virtual-circle center
    p_0dot = R_nb(1:2,1:2) * R_x2 * H * Lambda ...
        * x_tilde - k_0 * positionError;

    % Store the response before advancing to sample k+1
    wopcData.nu(i,:) = nu';
    wopcData.eta(i,:) = eta';
    wopcData.p0(i,:) = p_0';
    wopcData.positionError(i,:) = positionError';
    wopcData.xError(i,:) = x_tilde';
    wopcData.integralError(i,:) = e_I';
    wopcData.sigma(i,:) = sigma';
    wopcData.n(i,:) = n';
    wopcData.nCommand(i,:) = n_c';
    wopcData.tauControl(i,:) = tau';
    wopcData.tauDisturbance(i,:) = w_body';

    % Euler integration to sample k+1 (Fossen 2027, Appendix B)
    nu = nu + h * nu_dot;                           % Forward Euler
    eta = eta + h * R_nb * nu;                      % Backward Euler
    e_I = e_I + h * e_Idot;
    p_0 = p_0 + h * p_0dot;
end

% ------------------------------------------------------------------------------
%% PLOTS
% ------------------------------------------------------------------------------
plotWOPC(wopcData,R_d,p_d);

end

% ------------------------------------------------------------------------------
function draw_arrow(startpoint,endpoint,headsize)
% Draw a filled arrow between two [x y] plot coordinates
v1 = headsize * (startpoint-endpoint) / 2.5;
theta = 22.5*pi/180;
rotMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotMatrix1 = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)];
v2 = v1 * rotMatrix;
v3 = v1 * rotMatrix1;
x1 = endpoint;
x2 = x1 + v2;
x3 = x1 + v3;
fill([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[0 0 0])
plot([startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)], ...
    'LineWidth',2,'Color',[0 0 0])
end

% -------------------------------------------------------------------------
function plotWOPC(wopcData,R_d,p_d)
% plotWOPC plots the WOPC trajectory, heading, and virtual-circle center

t = wopcData.time / 60;           % Plot time in minutes
x = wopcData.eta(:,1);
y = wopcData.eta(:,2);
psi = wopcData.eta(:,3);
x0 = wopcData.p0(:,1);
y0 = wopcData.p0(:,2);

%% FIGURE 1: HEADING AND CIRCLE-CENTER COORDINATES
figure(1); clf
subplot(211)
plot(t,rad2deg(psi),'-.r','LineWidth',2)
hold on
% The bow points into the incoming weather, opposite its travel direction
psiOptimal = rad2deg(wopcData.headingEquilibrium);
plot([t(1),t(end)],[psiOptimal,psiOptimal],'-k','LineWidth',2)
hold off
title('Weather-optimal heading (deg)','FontSize',14)
xlabel('Time (min)','FontSize',12)
legend('Heading angle, \psi','Optimal heading','Location','best')
grid on

subplot(212)
plot(t,x0,'-.r',t,y0,'b-','LineWidth',2)
title('Virtual-circle center coordinates (m)','FontSize',14)
xlabel('Time (min)','FontSize',12)
legend('x_0 (m)','y_0 (m)','Location','best')
grid on

%% FIGURE 2: NORTH-EAST TRAJECTORY AND VIRTUAL CIRCLES
figure(2); clf
hold on

% Constant mean-load arrow in the North-East plot.
loadPlot = [wopcData.loadNED(2),wopcData.loadNED(1)];
loadPlot = loadPlot / norm(loadPlot);
arrowStart = [70,0];
arrowEnd = arrowStart + 15 * loadPlot;
draw_arrow(arrowStart,arrowEnd,0.5);
text(mean([arrowStart(1),arrowEnd(1)]), ...
     mean([arrowStart(2),arrowEnd(2)])+5,'w','fontsize',13)

% Ship trajectory and the initial virtual circle centered at p_d
plot(y(1:3:end),x(1:3:end),'k','LineWidth',1.5)
theta = 0:0.01:2*pi;
hVirtualCircle = plot(p_d(2) + R_d*sin(theta), ...
    p_d(1) + R_d*cos(theta), ...
    'r','LineWidth',1.5);

% Ship marker geometry.
L = 7;
B = 1.5;
YY = [0 B/2 B/2 -B/2 -B/2 0];
XX = [L/2 L/3.5 -L/2 -L/2 L/3.5 L/2];
Ao = atan2(YY,XX);
R = hypot(XX,YY);

% Select approximately equal arc-length positions instead of equal time
% intervals. This avoids clustering ship outlines after the vessel slows
nShipMarkers = 20;
pathLength = [0; cumsum(hypot(diff(x),diff(y)))];
if pathLength(end) > 0
    markerDistance = linspace(0,pathLength(end),nShipMarkers);
    shipIndex = zeros(1,nShipMarkers);
    for k = 1:nShipMarkers
        [~,shipIndex(k)] = min(abs(pathLength-markerDistance(k)));
    end
    shipIndex = unique(shipIndex,'stable');
else
    shipIndex = 1;
end

% Intermediate ships.
intermediateIndex = setdiff(shipIndex,[1,length(t)],'stable');
for i = intermediateIndex
    A = Ao + psi(i);
    YN = y(i) + R .* sin(A);
    XN = x(i) + R .* cos(A);
    plot(YN,XN,'Color','blue','LineWidth',1.5)
end

% Fill the initial and final ships so the endpoints stand out
endpointIndex = unique([1,length(t)]);
for i = endpointIndex
    A = Ao + psi(i);
    YN = y(i) + R .* sin(A);
    XN = x(i) + R .* cos(A);
    fill(YN,XN,'k','EdgeColor','k','LineWidth',2.0)
end

% Show six snapshots spaced along center travel, so the transient circle
% motion remains visible even when most of the simulation is at equilibrium.
centerDistance = [0; cumsum(hypot(diff(x0),diff(y0)))];
circleDistance = linspace(0,centerDistance(end),6);
circleIndex = zeros(1,6);
for k = 1:6
    [~,circleIndex(k)] = min(abs(centerDistance-circleDistance(k)));
end
circleIndex(end) = length(t);
circleIndex = unique(circleIndex,'stable');
% The initial circle is already drawn above.
for k = 2:length(circleIndex)
    i = circleIndex(k);
    plot(y0(i) + R_d*sin(theta),x0(i) + R_d*cos(theta), ...
        'r','LineWidth',1.0,'HandleVisibility','off')
end
hCenterPath = plot(y0,x0,'k--','LineWidth',1.5);

legend([hVirtualCircle,hCenterPath], ...
    {'Virtual-circle snapshots','Circle-center trajectory'},'FontSize',11,'Location','northeast')

hold off
xlabel('East (m)','FontSize',12)
ylabel('North (m)','FontSize',12)
title('Weather-optimal positioning response','FontSize',16)
grid on
axis equal
axis([-80 140 -80 80]);

end
