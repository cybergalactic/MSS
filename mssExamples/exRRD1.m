% exRRD1: Roll and sway-yaw transfer functions for the Son and Nomoto (1982)
% container ship.
%
% References:  
%   K. Son og K. Nomoto (1982). On the Coupled Motion of Steering and 
%   Rolling of a High Speed Container Ship, Naval Architect of 
%   Ocean Engineering 20:73-83. From J.S.N.A., Japan, Vol. 150, 1981.
% 
% Author:    Thor I. Fossen
% Date:      2001-10-31
% Revisions: 

U=7.0;    % Service speed

% Normalization variables
rho = 1025;                 % Water density (kg/m^3)
L = 175;                    % Length of ship (m)
 
% Lineari model using nondimensional matrices and states with dimension; 
% see Eq. (D.13) in Fossen (2021, Appendix D): 
%   T * M' * inv(T) * v_dot + (U/L) * T * N' * inv(T) * v + 
%      (U/L)^2 * T * G' * inv(T) * eta = (U^2/L) * T * b' * delta
%   where eta = [x y psi]' and nu = [v p r]'
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

% State-space model
Minv = inv(T*M*Tinv);
A11 = - Minv * (U/L)*T*N*Tinv;
A12 = - Minv * (U/L)^2*T*G*Tinv;
B1  =   Minv * (U^2/L)*T*b;

A = [ A11        A12(:,2:3)
      0   1   0   0   0
      0   0   1   0   0    ]
  
B = [ B1 ; 0 ; 0 ]

% Roll transfer function (removing yaw integrator when using model reduction)
roll = ss(A(1:4,1:4),B(1:4,1),[0 0 0 1],0)
yaw  = ss(A(1:4,1:4),B(1:4,1),[0 0 1 0],0)
yaw_integrator = tf(1,[1 0]);

% Decoupled reduced-order models
red_yaw  = ss(modred(yaw,[2,4],'del'));
red_roll = ss(modred(roll,3,'del'));

w = logspace(-3,0);
[mag1,phase1] = bode(series(yaw,yaw_integrator),w);
[mag2,phase2] = bode(series(red_yaw,yaw_integrator),w);
[mag3,phase3] = bode(roll,w);
[mag4,phase4] = bode(red_roll,w);

figure(gcf)
subplot(211),semilogx(w,20*log10(mag1(:)),':'),grid
xlabel('Frequency [rad/s]'),title('Gain [dB]')
hold on
semilogx(w,20*log10(mag2(:)))
semilogx(w,20*log10(mag3(:)),'linewidth',2)
semilogx(w,20*log10(mag4(:)))
hold off
axis([0.001 1 -40 10]);
legend('Yaw model','Decoupled yaw model','Roll model','Decoupled roll model')

subplot(212),semilogx(w,phase1(:),':'),grid
xlabel('Frequency [rad/s]'),title('Phase [deg]')
hold on
semilogx(w,phase2(:))
semilogx(w,phase3(:)-180,'linewidth',2)
semilogx(w,phase4(:))
hold off

set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

% Display transfer functions
zpk(series(yaw,yaw_integrator))
zpk(roll)

% Display reduced-order transfer functions
zpk(series(red_yaw,yaw_integrator))
zpk(red_roll)
