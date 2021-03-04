% Iserwood72 
%
% Data file containing the wind coefficients CX, CY, and CN for merchant ships. 
% Generated using the formulas of Isherwood (1972). 
%
% Variables:
% gamma_r = relative wind angle (rad)
% V_r     = relative wind speed (m/s)
% L	      = length overall (m)
% B	      = beam (m)
% A_L     = lateral projected area (m^2)
% A_T     = transverse projected area (m^2)
% A_SS    = lateral projected area of superstructure (m^2)
% S	      = length of curvature surrounding thelateral projection of model
% 		    minus the length of the waterline. Also exclude slender bodies such as masts and ventilators (m)
% C	      = positive distance from bow to the centroid of the lateral projected area (m)
% M	      = number of distinct groups of masts or king posts seen in lateral projection; king posts close 
%           against the bridge front are not included
%
% Author:    Thor I. Fossen, Marine Cybernetics AS
% Date:      18th March 2003
% Revisions: 
%
% Wind coeffisients / implemented in Simulink using interpolation in the CX,CY
% and CN tables
%
% cx =   A0 + A1*2*A_L/L^2 + A2*2*A_T/B^2 + A3*(L/B) + A4*(S/L) + A5*(C/L) + A6*M;
% cy = -(B0 + B1*2*A_L/L^2 + B2*2*A_T/B^2 + B3*(L/B) + B4*(S/L) + B5*(C/L) + B6*A_SS/A_L);
% cn = -(C0 + C1*2*A_L/L^2 + C2*2*A_T/B^2 + C3*(L/B) + C4*(S/L) + C5*(C/L));
% 
% Wind forces and moment
%
% tauX = 0.5*cx*rho_a*V_r^2*A_T;
% tauY = 0.5*cy*rho_a*V_r^2*A_L;
% tauN = 0.5*cn*rho_a*V_r^2*A_L*L;
% 
% tau_w = [tauX,tauY,tauN]';
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

% constant
rho_a = 1.224;              % density of air at 20 C

% relative wind directions for CX, CY and CN data
gamma_r_iserwood = (pi/180)*(0:10:180)';   % rad

% CX_data = [A0	A1	A2	A3	A4	A5	A6]
CX_data= [...  
 	2.152	-5.00	0.243	-0.164	0	    0       0	
 	1.714	-3.33	0.145	-0.121	0 	    0	    0	
 	1.818	-3.97	0.211	-0.143	0	    0	    0.033	
 	1.965	-4.81	0.243	-0.154	0	    0   	0.041	
 	2.333	-5.99	0.247	-0.190	0    	0	    0.042
 	1.726	-6.54	0.189	-0.173	0.348	0	    0.048	
 	0.913	-4.68	0	    -0.104	0.482	0	    0.052	
 	0.457	-2.88	0	    -0.068	0.346	0	    0.043	
 	0.341	-0.91	0	    -0.031	0	    0    	0.032	
 	0.355	0	    0	    0   	-0.247	0	    0.018	
 	0.601	0	    0	    0	    -0.372	0	    -0.020
 	0.651	1.29	0	    0	    -0.582	0	    -0.031	
 	0.564	2.54	0	    0	    -0.748	0	    -0.024	
 	-0.142	3.58	0	    0.047	-0.700	0	    -0.028	
 	-0.677	3.64	0	    0.069	-0.529	0	    -0.032	
 	-0.723	3.14	0	    0.064	-0.475	0	    -0.032	
	-2.148	2.56	0	    0.081	0	    1.27	-0.027	
	-2.707	3.97	-0.175	0.126	0	    1.81	0	
 	-2.529	3.76	-0.174	0.128	0    	1.55	0	      ];

% CY_data = [B0	B1	B2	B3	B4	B5	B6]
CY_data = [...
    0       0       0       0       0       0       0       
 	0.096	0.22	0	    0	    0   	0       0   	
 	0.176	0.71	0	    0	    0	    0	    0   	
 	0.225	1.38	0	    0.023	0	    -0.29	0   	
 	0.329	1.82	0	    0.043   0   	-0.59	0   	
 	1.164	1.26	0.121	0	    -0.242	-0.95	0   	
 	1.163	0.96	0.101	0	    -0.177	-0.88	0   	
 	0.916	0.53	0.069	0	    0   	-0.65	0   	
 	0.844	0.55	0.082	0	    0   	-0.54	0   	
 	0.889	0	    0.138	0	    0   	-0.66	0   
 	0.799	0	    0.155	0	    0    	-0.55	0   	
 	0.797	0	    0.151	0	    0	    -0.55	0   	
 	0.996	0	    0.184	0	    -0.212	-0.66	0.34	
 	1.014	0	    0.191	0	    -0.280	-0.69	0.44	
 	0.784	0	    0.166	0	    -0.209	-0.53	0.38	
 	0.536	0	    0.176	-0.029	-0.163	0	    0.27	
 	0.251	0	    0.106	-0.022	0	    0	    0	    
 	0.125	0	    0.046	-0.012	0	    0	    0	    
    0       0       0       0       0       0       0     ];

% CN_data = [C0	C1	C2	C3	C4	C5]
CN_data = [...
    0       0       0       0       0       0       
 	0.0596	0.061	0	    0	    0	    -0.074	
 	0.1106	0.204	0	    0	    0	    -0.170	
 	0.2258	0.245	0   	0	    0	    -0.380	
 	0.2017	0.457	0	    0.0067	0	    -0.472	
 	0.1759	0.573	0	    0.0118	0	    -0.523	
 	0.1925	0.480	0	    0.0115	0	    -0.546	
 	0.2133	0.315	0	    0.0081	0	     -0.526	
 	0.1827	0.254	0	    0.0053	0	    -0.443	
 	0.2627	0	    0	    0	    0	    -0.508	
 	0.2102  0	    -0.0195	0	    0.0335	-0.492	
 	0.1567	0	    -0.0258	0	    0.0497	-0.457
 	0.0801	0	    -0.0311	0	    0.0740	-0.396	
 	-0.0189	0	    -0.0488	0.0101	0.1128	-0.420	
 	0.0256	0	    -0.0422	0.0100	0.0889	-0.463
 	0.0552	0	    -0.0381	0.0109	0.0689	-0.476
 	0.0881	0	    -0.0306	0.0091	0.0366	-0.415
 	0.0851	0	    -0.0122	0.0025	0   	-0.220	
    0       0       0       0       0       0       ];
