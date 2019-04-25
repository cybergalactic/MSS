% viscous compensation
B44v = 1e10

% constants
rho = 1025
g   = 9.81

% WAMIT
%  C(3,3),C(3,4),C(3,5):   973.14       0.0000       0.0000    
%  C(4,4),C(4,5),C(4,6):               0.26303E+06   0.0000       0.0000    
%         C(5,5),C(5,6):                            0.26303E+06   0.0000 
C33 =  973.14*rho*g
C44 = 0.20429e6*rho*g
C55 = 0.20429e6*rho*g

% main data
m   =  (50708.6   +   50713.5   +   50713.9 )*rho/3

% radii of gyration
GM  = 4.0
R44 = 33.0     
R55 = 34.0
R66 = 37.5

% mass properties
xg = 0
yg = 0
zg = -0.5;
r_g  = [xg yg zg]
I_CO =  m*diag([R44^2 R55^2 R66^2 ])

% mass matrix about body axes (CO)
M_RB = [ ...
    m*eye(3)       -m*Smtrx(r_g)
    m*Smtrx(r_g)   I_CO]

% CO DATA WAMIT
T_p = 10
w_p = 2*pi/T_p
Ix  = I_CO(1,1)

%     ADDED-MASS AND DAMPING COEFFICIENTS
%      I     J         A(I,J)         B(I,J)
% 
%      1     1   1.344507E+04   1.537923E+01
%      1     5  -2.142362E+05  -4.468336E+01
%      2     2   4.644635E+04   3.484745E+01
%      2     4   6.194999E+05   2.163600E+02
%      3     3   6.497232E+04   6.010053E+02
%      4     2   6.198119E+05   2.160812E+02
%      4     4   6.960710E+07   1.341632E+03
%      5     1  -2.152238E+05  -4.458025E+01
%      5     5   5.548998E+07   1.295613E+02
%      6     6   5.011212E+07   5.294080E+00
A11 = 1.308480e4*rho;
A22 = 4.564919e4*rho;      
A66 = 4.975603e7*rho;

% 20 and 40 s data
A33 = 5.961171e4*rho;
A44 = 6.960710e7*rho;
A55 = 6.107354e7*rho;

B44  = 1.341632e3*1025*w_p

ratio = [A11/m A22/m A33/m A44/I_CO(1,1) A55/I_CO(2,2) A66/I_CO(1,1)]
    
zeta_p = B44/(2*(Ix+A44)*w_p)
zeta_v = B44v/(2*(Ix+A44)*w_p)

zeta_total = zeta_p + zeta_v

% transform added mass to CO
M_A = diag([A11 A22 A33 A44 A55 A66]);

% total mass in CO
M  = M_RB + M_A;

% DP system
T_surge        = 125;   % supply 100-120, tanker 150-200
T_sway         = 155;
T_yaw          = 91;
zeta_PD        = 0.7;          % Kongsberg setting

w_surge = 2*pi/T_surge;
w_sway  = 2*pi/T_sway;
w_yaw   = 2*pi/T_yaw;

GM_T    = C44/(m*g)
T_heave = 2*pi*sqrt(M(3,3)/C33)
T_roll  = 2*pi*sqrt(M(4,4)/C44)
T_pitch = 2*pi*sqrt(M(5,5)/C55)


% wamit axes
MRB_CO = diag([m m m m*R44^2 m*R55^2 m*R66^2]);
MRB_CO(1,5) = m*zg;
MRB_CO(5,1) = MRB_CO(1,5);
MRB_CO(2,4) = -m*zg;
MRB_CO(4,2) = MRB_CO(2,4);
MRB_CO(3,5) = -m*xg;
MRB_CO(5,3) = MRB_CO(3,5);

format long g
MRB = sprintf('%1.2f %1.2f %1.2f %1.2f %1.2f %6.2f\n',MRB_CO)


