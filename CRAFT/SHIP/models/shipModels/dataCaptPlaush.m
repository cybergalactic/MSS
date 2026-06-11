%  Ship's particulars data: Captain Plaushevskiy (General Cargo)
% 
% 
% Author:    Aleksandr D. Pipchenko
% Date:      10th April 2007
clc;

%  Dimensions
 L = 140;               %Length between perpendiculars, m
 B = 20.6;              %Breadth overall, m
 Ta = 9;                %Draught Aft, m
 Tf = 8.95;             %Draught Fore, m
 Tm = (Ta+Tf)/2;        %Midle Draught, m
 i = 17;                %Number of theoretial frame from U to V form  of the hull
 
 Ac = 80;               %Aft underhull area, m^2
 Hc = 0.6886;           %Hull coefficient
 x = 0.5;               %Hull influence coefficient
 %Rudder
 Ar = 17;               %Rudder Area, m^2
 Ar0 = 3.75;            %Rudder area not under propeller's jet, m^2
 Ardp = 13.25;          %Rudder area under propeller's jet, m^2
 hr = 6.8;              %Spindle Rudder Height, m
 Ir = 70;               %Distance from 10-th theoretical frame to rudder, m 
 cr = 1;                %Rudder shape correcting coefficient (for rectangular or trapeze rudders cr = 1)
 k = 1;                 %Addition coefficient (For common rudders k=1, for rudders after ruderpost k = 1.3, for tail rottors k = 0.9)
 %Propeler
 Dp = 5.3;              %Propellers diameter, m
 Pitch = 4.998;         %Propellers pitch, m
 P_Dp = (Pitch/Dp);     %Pitch Ratio
 Z = 4;                 %Number of Blades
 Q = 0.6;               %Disk Ratio (Expanded Area)
%Miscelanous
 Ctv = 0.005;           %coarse pitch stop coefficient
 
 %Water resistence forces coefficients
  
 Cbyb = 0.0921;         
 c2 = 0.6101;           
 c3 = 0.064;
 m1 = 0.0565;    
 m3 = 0.0206;
 m4 = 0.0134;
 
 %Constants
 rho = 1025;            %Density of sea water, kg/m^3
 g = 9.8;               %Gravity constant, N/kg
 CA = pi/180;           %Degrees to radians constant
 
  %1-st step formulas
 psi1 = (Ta-Tf)/L;      %Static trim
 psi2 = 0;              %Dynamic trim
  Cappa = 1 - 3*Ac/(L*Tm*(20-i))+0.054*(psi1+psi2)*L/Tm; %Corrected hull coefficient
  Ald = L*Tm*Cappa;     %Corrected area of deepped logitudinal plane's part, m^2 
  Are = 1.05*Ar;        %Effective rudder area, m^2
  
  Xdp=(Ar0+Ardp*sqrt(1+Ctv))/(Ar0+Ardp*(1+Ctv)); %Propellers influence coefficient
  Xe = 2*x*Xdp;   %effective coefficient of propeller's & hull influence on rudder

 Caryr = 2*pi*cr/(1+2*Ar/hr^2);
 
  Cayr = 4*k*Caryr;     %Transverse rudder force on attack angle coefficient
 
  Irm = Ir/L;           %Relative distance from 10-th theoretical frame to rudder
  
 m2 = 0.0271;           %Water resistance coefficient 
 Cwmw =(0.739+8.7*Tm/L)*(1.611*Cappa^2-2.873*Cappa+1.33); %Demphing moment coeffisient
 
 %Kinematic parameters
 
 V = B*Tm*Hc*L;         %Ships current volume displacement, m^3
 D = 18321000;          %Ship's current displacement, kg
 k11 = 0.08;
 k22 =0.042;	        %Connected masses coefficients 
 k66 =0.032;
 Izz = 0.05*D*L^2;      %Vertical axes's moment of inertia
 
 %Control Forces (Propeller)
 
 Kr = 1.5962;
 Kp = ((Q*Z)^(1/3))*(0.225*sin(P_Dp)^2+0.098*sin(P_Dp));%Thrust coefficient
 Cx0 = 0.0485;