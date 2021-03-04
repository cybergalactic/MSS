% ExRRD2     Rudder-Roll Damping (RRD) system for the Son and Nomoto container ship
%            Calls Son and Nomoto model with wave-induced roll disturbance: see RRDcontainer.m
% Author:    Thor I. Fossen
% Date:      31 October 2001
% Revisions: 

conversion % load conversion factors

Ns = 20000; % no. of samples
h  = 0.05;  % sample time

% Normalization variables
rho = 1025;                 % water density (kg/m^3)
L = 175;                    % length of ship (m)
U0 = 7.3;

% -----------------------------------------------------------------------------------------
% SHIP MODEL
% -----------------------------------------------------------------------------------------
 
% Linear model using nondimensional matrices and states with dimension (see Lcontainer.m): 
% TM'inv(T) dv/dt + (U/L) TN'inv(T) v + (U/L)^2 TG'inv(T) eta = (U^2/L) T b' delta

% nu = [v p r]
T    = diag([ 1 1/L 1/L]);
Tinv = diag([ 1 L L ]);

M = [ 0.01497    0.0003525     -0.0002205       
     -0.0002205  0.0000210      0              
      0.0003525  0              0.000875  ];

N = [ 0.012035   0.00522    0  
     -0.000314   0.0000075  0.0000692      
      0.0038436  -0.000213  0.00243    ];

G = [ 0 0.0000704  0
      0 0.0004966  0 
      0 0.0001468  0];

b = [-0.002578 0.0000855 0.00126  ]';

% ship state-space model
Minv = inv(T*M*Tinv);
A11 = - Minv * (U0/L)*T*N*Tinv;
A12 = - Minv * (U0/L)^2*T*G*Tinv;
B1  =   Minv * (U0^2/L)*T*b;

A = [ A11        A12(:,2:3)
      0   1   0   0   0
      0   0   1   0   0    ]
  
B = [ B1 ; 0 ; 0 ]

[PHI,DELTA] = c2d(A,B,h);

C = [0 1 0 0 0
     0 0 1 0 0
     0 0 0 1 0
     0 0 0 0 1 ];

% LQ controller gains
Q = diag([10000 1000 10 1 ]);
R = 5;
[G1,G2] = lqtracker(A,B,C,Q,R)

% eigenvalues, damping ratios, natural frequencies
damp(A)
damp(A+B*G1)

% -----------------------------------------------------------------------------------------
% SIMULATION OF AUTOPILOT AND RRD SYSTEMS
% -----------------------------------------------------------------------------------------
n_ref = 70; % desired rpm

x = [U0 0 0 0 0 0 0 0 0 n_ref 0 0]';  % x = [u v r x y psi p phi delta n xw1 xw2]'
xout = zeros(Ns+1,length(x));  

wf1 = 0.1;    % low-pass filter frequnecy
wf2 = 0.05;   % high-pass filter frequency

xyawf  = [0 0 0]';  % filter states
xrollf = [0 0]';
xw     = [0 0]';    % wave states
% --------------------------- MAIN LOOP --------------------------------------------

for i=1:Ns+1, 
    
    % state vectors for yaw and roll subsystems 
    xyaw  = [x(2) x(3) x(6)]'; 
    xroll = [x(7) x(8)]';
   
    % discrete time low-pass filter
    phif1   = exp(-h*wf1);
    xyawf  = phif1*xyawf + (1-phif1)*xyaw;

    % discrete time high-pass filter 
    phif2   = exp(-h*wf2);
    yrollf = xrollf + xroll;
    xrollf = phif2*xrollf + (phif2-1)*xroll;
    
    % x = [v p r phi psi ]'
    x_fb = [xyawf(1) yrollf(1) xyawf(2) yrollf(2) xyawf(3) ]';  % filtered states

    psi_d = 10*D2R;  % desired heading
    
    % control law:  u_ref = G1 * x + G2 * C * x_d
    if i<6000,
        u_ref = [G1(1) 0 G1(3) 0 G1(5)]*x_fb + G2*[0 0 0 psi_d]';   % course-changing, no roll damping
    elseif i>6000 & i<15000,
        u_ref = G1*x_fb + G2*[0 0 0 psi_d]';                        % course-keeping, roll damping
    else
        u_ref = [G1(1) 0 G1(3) 0 G1(5)]*x_fb + G2*[0 0 0 psi_d]';   % course-keeping, no roll damping    
    end 
   
    % save simulation results
    xout(i,:) = x';    

    % nonlinear ship model with wave disturbance
    [xdot,U] = RRDcontainer(x,[u_ref; n_ref]);
    x = x + h*xdot;  
  
end

% -----------------------------------------------------------------------------------------
% PLOTS
% -----------------------------------------------------------------------------------------

t = h*(0:1:Ns)';
figure(2)
subplot(311)
plot(t,R2D*xout(:,6)); hold on
plot(t,R2D*psi_d*ones(Ns+1,1),'linewidth',2); hold off; title('yaw angle \psi (deg)')
subplot(312),plot(t,R2D*xout(:,8)); title('roll angle \phi (deg)')
subplot(313),plot(t,R2D*xout(:,9)); title('rudder angle \delta (deg)')