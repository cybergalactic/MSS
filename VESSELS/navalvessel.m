function [xdot,U] = navalvessel(x, tau)
% Noninear maneuvering model in surge, sway, roll and yaw for a multipurpose
% naval vessel. The surge equation is decoupled except a centripetal term.
%
% [xdot,U] = navalvessel(x, tau)
% 
% Inputs:
% x = [u v p r phi psi]'
% tau = [Xe Ye Ke Ne]' 
%
% where
% u     = surge velocity          (m/s)
% v     = sway velocity           (m/s)
% p     = roll velocity           (rad/s)
% r     = yaw velocity            (rad/s)
% phi   = roll angle              (rad)
% psi   = yaw angle               (rad)
%
% Xe is the surge external force (e.g. rudder and fin forces)
% Ye is the sway external force  
% Ke is the sway external force  
% Ne is the sway external force   
%
% Reference: Blanke M. and Christensen A. (1993) "Rudder-roll 
% damping autopilot robustness to sway-yaw-roll couplings." 
% 10th Ship Control Systems Symposium, Ottawa, Canada.
%  
% Notes: 1 - The model does not include rudder Machinery.
%        2 - The parameters of the model should be defined 
%            in the structures h and const before using the function.  
%
% Author:     Tristan Perez
% Revisions:
%   - Original model from A. G. Jensen and M.S.Chislett (Danish Maritime
%     Institute) 1983-1989. 
%	- Adapted for Matlab by Mogens Blanke and Antonio Tiano 1996
%   - Modified for Matlab 5.3 implementation by Mogens Blanke 1997
%   - Modified for Simulink by Mogens Blanke and Tristan Perez, 2001 
%   - This version is further modified by Tristan Perez to match the data 
%     of the vessel design of ADI-Limited Australia.

% Vessel Data
const.rho_water     =	1014.0;	        %	water density	[kg/m^3]	
const.rho_air		=	1.225	;	    %	air density		[kg/m^3]	
const.g				=	9.81;	        %	gravity constant	[m/s^2]	
const.deg2rad 		=	pi/180;	        %	degrees to radians	
const.rad2deg 		=	180/pi;	        %	rad to degrees		
const.ms2kt			=	3600/1852;	    % 	m/s to kt 			
const.kt2ms 		=	1852/3600;	    %	kt to m/s			
const.RPM2rads		=	2*pi/60;	    %	RPM to rad/s		
const.rads2RPM		=   60/[2*pi];	    %	rad/s to RPM		
const.HP2W			=	745.700;	    %	HP to Watt			

% Struct rudder  (Modified by T.Perez)
rudder.sp    =1.5;                   % span
rudder.A     =1.5;                   % Area
rudder.ar    =3;                     % aspect ratio
rudder.dCL   =0.054; % 1/deg         % dCL/d a_e
rudder.stall =23;                    % a_stall 

% Main Particulars (Modified by T.Perez)
h.Lpp    =  51.5 ;                  % Length between perpendiculars [m]
h.B      =  8.6  ;                  % Beam over all  [m]
h.D	     =  2.3  ;                  % Draught [m]     

%Load condition (Modified by T.Perez)
h.disp   =  357.0;                   % Displacement  [m^3]
h.m      =  h.disp*const.rho_water;  % Mass [Kg]
h.Izz    =  47.934*10^6 ;            % Yaw Inertia
h.Ixx    =  2.3763*10^6 ;            % Roll Inertia
h.U_nom  =  8.0   ;	                 % Speed nominal [m/sec] (app 15kts) 
h.KM		=  4.47;	             %  [m] Transverse metacentre above keel
h.KB		=  1.53;	             %  [m] Transverse centre of bouancy
h.gm 		=  1.1;	                 %  [m]	Transverse Metacenter
h.bm 		=  h.KM - h.KB;
h.LCG       = 20.41 ;                % [m]	Longitudinal CG (from AP considered at the rudder stock)
h.VCG       = 3.36  ;                % [m]	Vertical  CG  above baseline
h.xG        = -3.38  ;               % coordinate of CG from the body fixed frame adopted for the PMM test  
h.zG  	    = -(h.VCG-h.D);          % coordinate of CG from the body fixed frame adopted for the PMM test  
h.m_xg	    = h.m * h.xG;
h.m_zg	    = h.m * h.zG;
h.Dp        = 1.6 ;                   % Propeller diameter [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The hydrodynamic derivatives are given in dimensional form, and follow
% from the original publication of Blanke and Christensen 1993.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for surge equation 
h.Xudot  	= -17400.0 ;
h.Xuau     	= -1.96e+003 ;
h.Xvr    	=  0.33 * h.m ;
   
% Hydrodynamic coefficients in sway equation
h.Yvdot = -393000 ; 
h.Ypdot = -296000 ; 
h.Yrdot = -1400000 ; 
h.Yauv  = -11800 ; 
h.Yur   =  131000 ; 
h.Yvav  = -3700 ; 
h.Yrar   =  0 ;
h.Yvar  = -794000 ; 
h.Yrav  = -182000 ; 
h.Ybauv =  10800 ; % Y_{\phi |v u|}
h.Ybaur =  251000 ; 
h.Ybuu  = -74 ; 


% Hydrodynamic coefficients in roll equation
h.Kvdot =  296000 ;
h.Kpdot = -774000 ;
h.Krdot =  0 ;
h.Kauv  =  9260 ;
h.Kur   = -102000 ;
h.Kvav  =  29300 ;
h.Krar  =  0 ;
h.Kvar  =  621000 ;
h.Krav  =  142000 ;
h.Kbauv =  -8400 ;
h.Kbaur =  -196000 ;
h.Kbuu  =  -1180 ;
h.Kaup  =  -15500 ;
h.Kpap  =  -416000 ;
h.Kp    =  -500000 ;
h.Kb    =  0.776*h.m*const.g;
h.Kbbb  =  -0.325*h.m*const.g ;

% Hydrodynamic coefficients in yaw equation
h.Nvdot =  538000 ;
h.Npdot =  0 ;
h.Nrdot = -38.7e6;
h.Nauv  = -92000 ;
h.Naur  = -4710000 ;
h.Nvav  =  0 ;
h.Nrar  = -202000000 ;
h.Nvar  =  0 ;
h.Nrav  = -15600000 ;
h.Nbauv = -214000 ;
h.Nbuar = -4980000 ;
h.Nbuau = -8000 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rename inputs of the function
u   = x(1);
v   = x(2);	 
p  	= x(3);   	
r  	= x(4); 
b  	= x(5); % b -denoted phi roll 		
psi = x(6);
U = sqrt(u^2+v^2);

% External forces
Xe  = tau(1);
Ye  = tau(2);
Ke  = tau(3);
Ne  = tau(4);

% Auxiliary variables
au = abs(u);
av = abs(v); 
ar = abs(r); 
ap = abs(p); 
ab = abs(b);
L2 = h.Lpp^2; 

% Total Mass Matrix 
M =[ (h.m-h.Xudot)  0   0   0   0   0;
   0 (h.m-h.Yvdot) -(h.m*h.zG+h.Ypdot) (h.m*h.xG-h.Yrdot) 0 0;
   0 -(h.m*h.zG+h.Kvdot) (h.Ixx-h.Kpdot) -h.Krdot 0 0;
   0 (h.m*h.xG-h.Nvdot) -h.Npdot (h.Izz-h.Nrdot) 0 0;
   0 0 0 0 1 0; 
   0 0 0 0 0 1] ;
% Hydrodynamic forces without added mass terms (considered in the M matrix)
Xh  = h.Xuau*u*au+h.Xvr*v*r;

Yh = h.Yauv*au*v + h.Yur*u*r + h.Yvav*v*av + h.Yvar*v*ar + h.Yrav*r*av ...
   + h.Ybauv*b*abs(u*v) + h.Ybaur*b*abs(u*r) + h.Ybuu*b*u^2;

Kh = h.Kauv*au*v +h.Kur*u*r + h.Kvav*v*av + h.Kvar*v*ar + h.Krav*r*av ...
   + h.Kbauv*b*abs(u*v) + h.Kbaur*b*abs(u*r) + h.Kbuu*b*u^2 + h.Kaup*au*p...
   + h.Kpap*p*ap +h.Kp*p +h.Kbbb*b^3-(const.rho_water*const.g*h.gm*h.disp)*b;

Nh = h.Nauv*au*v + h.Naur*au*r + h.Nrar*r*ar + h.Nrav*r*av...
   +h.Nbauv*b*abs(u*b) + h.Nbuar*b*u*ar + h.Nbuau*b*u*au;
 
% Rigid-body centripetal accelerations
Xc =   h.m*(r*v+h.xG*r^2-h.zG*p*r);  
Yc = - h.m*u*r;
Kc =   h.m*h.zG*u*r;
Nc = - h.m*h.xG*u*r;

% Total forces
F1 = Xh+Xc+Xe;
F2 = Yh+Yc+Ye;
F4 = Kh+Kc+Ke;
F6 = Nh+Nc+Ne;
 
xdot = M\[F1; F2; F4; F6; p; r];




