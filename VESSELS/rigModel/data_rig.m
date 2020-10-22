% DATA_RIG    Configuration program for computation of LF semi-sub. model 
%             Edit file by inlcuding main semi-sub. dimensions from drawing.
%             All coordinates are given in COO by using SI-units.
%             
%             Inputs: main dimensions 
%             Ouput:  DATCONF   (array struct) used by rig.m
%
%
% _________________________________________________________________
%
%         L11      L12      L13          
%      |-------------------------|   Pi = pontoon i  (x)
%   P1 |     O     x o        o  |   Lij = leg j pontoon i (o)
%      |-------------------------|              
%                    |                 
%                    |       
%     ---------------------------------------------> x
%                    |
%                    |
%      |-------------------------|     
%   P2 |     o     x o        o  |    
%      |-------------------------|
%          L21     L22       L23    
%                    |
%                    |  
%                    v   y
%
% Author:     12 March 2000  Thor I. Fossen
% ________________________________________________________________

draft = 23.6;    

% Constants 
DATCONF.N_p  = 2;        % number of pontoons
DATCONF.N_l  = 4;        % number of legs

DATCONF.m    = 27162500;  % mass of vessel
DATCONF.rho  = 1026;     % sea water density
DATCONF.g    = 9.81;     % acceleration of gravity

DATCONF.r_x = 30.0;      % radius of gyration        
DATCONF.r_y = 32.0;
DATCONF.r_z = 37.0;

DATCONF.R_CG = [-0.5 0 -19.5]'; % center of gravity in COO

DATCONF.T_x = 100.0;      % time constants in surge, sway and yaw
DATCONF.T_y = 200.0;
DATCONF.T_n =  80.0;

%%%%%%%%%%%%%%%%%%%% OPTIONALLY PARAMETERS %%%%%%%%%%%%%%

% NOT IMPLEMENTED IN THIS VERSION - COMPUTED BY CONF_RIG

% set values to zero if unknown 

DATCONF.GMTtran   = 0;    % Transverse metasentric height, transit [m]
DATCONF.GMLtran   = 0;    % Longitudinal metacentric height, transit [m]

DATCONF.GMTball   = 0;    % Transverse metasentric height, ballasted (survival) [m]
DATCONF.GMLball   = 0;    % Longitudinal metacentric height, ballasted (survival) [m]

DATCONF.GMTload   = 0.5;    % Transverse metasentric height, fully loaded (operational) [m]
DATCONF.GMLload   = 5.4;    % Longitudinal metacentric height, ballasted (operational) [m]


%%%%%%%%%%%%%%%%%%%%% PONTOONS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Main dimensions
DATCONF.H_p = 13.0;                  % pontoon height
DATCONF.B_p = 12.0;                  % pontoon width
DATCONF.L_p = 84.6;                  % pontoon length

% pontoon block coefficients: cB = nabla_p/(Lp Bp Hp) <= 1.0
DATCONF.CB_p = 0.75;                 % volume correction for rounded corners
                           
% volume centers in COO (only xy-coordinates)
DATCONF.R_p(1,:) = [0  21.0 ];            % port pontoon
DATCONF.R_p(2,:) = [0 -21.0 ];            % starboard pontoon

%%%%%%%%%%%%%%%%%%%%% LEGS %%%%%%%%%%%%%%%%%%%%%%%%%%

% shape of legs
DATCONF.legs = 'box';              % legs = {'box' or 'cylinder'} 

% block coefficient for legs
DATCONF.CB_l = 0.80;               % for cylinder shaped legs: CBl = 1.0

% Main dimensions of legs
DATCONF.B_l(1) = 16.0;             % width  of leg (for cylinder B_l = diameter)
DATCONF.L_l(1) = 12.0;             % length of leg (for cylinder L_l = diameter)

DATCONF.B_l(2) = 16.0;
DATCONF.L_l(2) = 12.0;

DATCONF.B_l(3) = 16.0;
DATCONF.L_l(3) = 12.0;

DATCONF.B_l(4) = 16.0;
DATCONF.L_l(4) = 12.0;

% volume centers in COO (only xy-coordinates)
DATCONF.R_l(1,:) = [ 24.0  21.0 ];
DATCONF.R_l(2,:) = [ 24.0 -21.0 ];
DATCONF.R_l(3,:) = [-24.0 -21.0 ];
DATCONF.R_l(4,:) = [-24.0  21.0 ];

%%%%%%%%%%%%%%%%%%%%% BRACINGS %%%%%%%%%%%%%%%%%%%%%%%%%%

DATCONF.D_b = 1.7;          % average bracing diameter
DATCONF.L_b = 2*30;         % total length of bracing
DATCONF.R_b = [0 0 -12];    % (x,y,z) coordinates for volume center 

