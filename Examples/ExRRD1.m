% ExRRD2     Roll and sway-yaw transfer functions for the Son and Nomoto container ship
% Author:    Thor I. Fossen
% Date:      21 October 2001
% Revisions: 
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org

U=7.0;    % service speed

% Normalization variables
rho = 1025;                 % water density (kg/m^3)
L = 175;                    % length of ship (m)
 
% Linear model using nondimensional matrices and states with dimension (see Fossen 2002): 
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

%  state-space model
Minv = inv(T*M*Tinv);
A11 = - Minv * (U/L)*T*N*Tinv;
A12 = - Minv * (U/L)^2*T*G*Tinv;
B1  =   Minv * (U^2/L)*T*b;

A = [ A11        A12(:,2:3)
      0   1   0   0   0
      0   0   1   0   0    ]
  
B = [ B1 ; 0 ; 0 ]

% roll transfer function (removing yaw integrator when using model reduction)
roll=ss(A(1:4,1:4),B(1:4,1),[0 0 0 1],0)
yaw      = ss(A(1:4,1:4),B(1:4,1),[0 0 1 0],0)
yaw_integrator = tf(1,[1 0]);

% decoupled reduced order models
red_yaw  = ss(modred(yaw,[2,4],'del'));
red_roll = ss(modred(roll,[3],'del'));

w = logspace(-3,0);
[mag1,phase1] = bode(series(yaw,yaw_integrator),w);
[mag2,phase2] = bode(series(red_yaw,yaw_integrator),w);
[mag3,phase3] = bode(roll,w);
[mag4,phase4] = bode(red_roll,w);

subplot(211),semilogx(w,20*log10(mag1(:)),':b'),grid
xlabel('Frequency [rad/s]'),title('Gain [dB]')
hold on
semilogx(w,20*log10(mag2(:)),'b')
semilogx(w,20*log10(mag3(:)),'k','linewidth',2)
semilogx(w,20*log10(mag4(:)),'k')
hold off
axis([0.001 1 -40 10]);
legend('yaw model','decoupled yaw model','roll model','decoupled roll model')

subplot(212),semilogx(w,phase1(:),':b'),grid
xlabel('Frequency [rad/s]'),title('Phase [deg]')
hold on
semilogx(w,phase2(:),'b')
semilogx(w,phase3(:)-180,'k','linewidth',2)
semilogx(w,phase4(:),'k')
hold off

% transfer functions
zpk(series(yaw,yaw_integrator))
zpk(roll)

% decoupled transfer functions
zpk(series(red_yaw,yaw_integrator))
zpk(red_roll)
