% ExObsCtr   Observability and Controllability of Ships (see supply.m)
% Author:    Thor I. Fossen
% Date:      11 July 2002
% Revisions: 20 Dec 2009   minor modifications

format compact
disp('Controllabilty and observability of offShore supply vessel length 76 m')

% Normalization variables

L   =  76.2;           % length of ship (m)
g   =  9.8;            % acceleration of gravity (m/s^2)
m   = 6000e3;          % mass (kg)

T    = diag([1 1 L]);
Tinv = diag([1 1 1/L]);

% Model matricses
Mbis = [1.1274         0          0
             0    1.8902    -0.0744
             0   -0.0744     0.1278];

Dbis = [0.0358        0        0
             0        0.1183  -0.0124
             0       -0.0041   0.0308];
   
Athr = diag([-1/100 -1/100 -1/100]);

M = m*Tinv^2*(T*Mbis*Tinv)
D = m*sqrt(g/L)*Tinv^2*(T*Dbis*Tinv)

% state space model
A = [ zeros(3,3) eye(3)     zeros(3,3)
      zeros(3,3) -inv(M)*D  inv(M)
      zeros(3,3) zeros(3,3) Athr      ];

B = [zeros(3,3); zeros(3,3); -Athr ];

F = [ zeros(3,3) zeros(3,3) eye(3)
      zeros(3,3) zeros(3,3) zeros(3,3)  
      zeros(3,3) inv(M) -inv(M)*D ];

H = [ eye(3) zeros(3,3) zeros(3,3)];


n=rank(ctrb(A,B))

n=rank(obsv(F,H))

% augmented state space modelwith integral action
C = [ eye(3) zeros(3,3) zeros(3,3)];
Aa = [zeros(3,3) C
      zeros(9,3) A ];

Ba = [zeros(3,3)
      B ];

n=rank(ctrb(Aa,Ba))
