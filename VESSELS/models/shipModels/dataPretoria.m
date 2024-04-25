%  Ship's particulars data: Hanjin Pretoria (Container Carrier)
% 
% 
% Author:    Aleksandr D. Pipchenko
% Date:      10th April 2007
clc;

%  Dimensions
 L = 266;               %Length between perpendiculars, m
 B = 32.2;              %Breadth overall, m
 Ta = 8;                %Draught Aft, m
 Tf = 9.13;             %Draught Fore, m
 Tm = (Ta+Tf)/2;        %Middle Draught, m
 i = 17;                %Number of theoretial frame from U to V form of the hull
 Ac = 17*12;            %Aft underhull area, m^2
 Hc = 0.817;            %Hull coefficient
 Betta = 0.997;         %Midel frame coefficient
 Phi = Hc/Betta;        %Logitudinal hull coefficient
 x = 0.5;               %Hull influence coefficient
 %Rudder
 Ar = 53.6;             %Rudder Area, m^2
 Ar0 = 10;              %Rudder area not under propeller's jet, m^2
 Ardp = 43.6;           %Rudder area under propeller's jet, m^2
 hr = 10;               %Spindle Rudder Height, m
 Ir = 132.6;            %Distance from 10-th theoretical frame to rudder, m 
 cr = 1;                %Rudder's shape correcting coefficient (for rectangular or trapeze rudders cr = 1)
 k = 1;                 %Additional coefficient (For common rudders k=1, for rudders after ruderpost k = 1.3, for tail rottors k = 0.9)
 %Propeler
 Dp = 8.2;              %Propeller diameter, m
 Pitch = 8.331;         %Propeller pitch, m
 P_Dp = (Pitch/Dp);     %Pitch Ratio
 Z = 5;                 %Number of Blades
 Q = 0.6;               %Disk Ratio (Expanded Area)
 %Miscelanous
 Ctv = 0.005;           %coarse pitch stop coefficient
 
 %Water resistence forces coefficients
 Cbyb = 0.14;           
 c2 = 0.63;              
 c3 = 0.058;
 m1 = 0.055;            
 m3 = 0.022;
 m4 = -0.005;
 
 %Constants
 rho = 1025;            %Density of sea water, kg/m^3
 g = 9.8;               %Gravity constant, N/kg
 CA = 3.141592654/180;  %Degrees to radians constant
 
%1-st step formulas
 psi1 = (Ta-Tf)/L;      %Static trim
 psi2 = 0;              %Dynamic trim
  Cappa = 1 - 3*Ac/(L*Tm*(20-i))+0.054*(psi1+psi2)*L/Tm; %Corrected hull coefficient
  Ald = L*Tm*Cappa;     %Corrected area of deepped logitudinal plane's part, m^2 
  Are = 0.94*Ar;        %Are = Ar0 + Ardp*(1+Ctv) %Effective rudder area, m^2
  
  Xdp=(Ar0+Ardp*sqrt(1+Ctv))/(Ar0+Ardp*(1+Ctv)); %Propellers influence coefficient
  Xe = 2*x*Xdp;         %Effective coefficient of propeller's & hull influence on rudder

 Caryr = 2*pi*cr/(1+2*Ar/hr^2);
 
  Cayr = 4*k*Caryr;     %Transverse rudder force on attack angle coefficient
 
  Irm = Ir/L;           %Relative distance from 10-th theoretical frame to rudder
  
  m2 = -0.02;           %m2 = -(log(1.023*Cappa))/(11.6*Cappa-9.29); %Water resistance coefficient 
 Cwmw =(0.739+8.7*Tm/L)*(1.611*Cappa^2-2.873*Cappa+1.33); %Demphing moment coeffisient
 
 %Kinematic parameters
 
 V = B*Tm*Hc*L;         %Ships current volume displacement, m^3
 D = V*rho;             %Ship's current displacement, kg
 k11 = 0.02;            %(0.465*B/L-0.028)*2*Tm/B;%
 k22 = 0.042;           %(1-0.5*B/L)*2*Tm/B          %Connected masses coefficients
 k66 = 0.032;           %(1.03-1.76*B/L)*2*Tm/B
 Izz = 0.05*rho*B*Tm*Hc*L^3;         %Vertical axes's moment of inertia
 
 %Control Forces (Propeller)
 
 Kr = 1.5962;
 Kp = ((Q*Z)^(1/3))*(0.225*sin(P_Dp)^2+0.098*sin(P_Dp));%Thrust coefficient
 Cx0 = 0.0485;