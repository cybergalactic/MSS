function [tau_w,CX,CY,CN] = isherwood72(gamma_r,V_r,Loa,B,ALw,AFw,A_SS,S,C,M)
% [tau_w,CX,CY,CN] = isherwood72(gamma_r,V_r,Loa,B,ALw,AFw,A_SS,S,C,M) returns the the wind 
% force/moment vector w_wind = [tauX,tauY,tauN] and the optionally wind coeffisients
% cx,cy and cn for merchant ships using the formulas of Isherwood (1972). 
%
% INPUTS:
% gamma_r = relative wind angle (rad)
% V_r     = relative wind speed (m/s)
% Loa     = length overall (m)
% B	      = beam (m)
% ALw     = lateral projected area (m^2)
% AFw     = frontal projected area (m^2)
% A_SS    = lateral projected area of superstructure (m^2)
% S	      = length of perimeter of lateral projection of model (m)
% 		    excluding waterline and slender bodies such as masts and ventilators (m)
% C	      = distance from bow of centroid of lateral projected area (m)
% M	      = number of distinct groups of masts or king posts seen in lateral
% 		    projection; king posts close against the bridge front are not included
%
% Author:    Thor I. Fossen
% Date:      10th September 2001
% Revisions: 19.04.2004, changed velocity from knots to m/s. This was a bug
%            20.11.2008, changed name from windcoef to isherwood72, updated
%                        signs and notation to comply with Blendermann (1994).

if nargin~=10, error('the number of inputs must be 10');end

% conversions and constants
rho_a = 1.224;             % density of air at 20 C
gamma_r = gamma_r*180/pi;  % rad2deg

% CX_data = [gamma_r 	A0	A1	A2	A3	A4	A5	A6	]
CX_data= [...  
0	2.152	-5.00	0.243	-0.164	0	    0       0	
10	1.714	-3.33	0.145	-0.121	0 	    0	    0	
20	1.818	-3.97	0.211	-0.143	0	    0	    0.033	
30	1.965	-4.81	0.243	-0.154	0	    0   	0.041	
40	2.333	-5.99	0.247	-0.190	0    	0	    0.042
50	1.726	-6.54	0.189	-0.173	0.348	0	    0.048	
60	0.913	-4.68	0	    -0.104	0.482	0	    0.052	
70	0.457	-2.88	0	    -0.068	0.346	0	    0.043	
80	0.341	-0.91	0	    -0.031	0	    0    	0.032	
90	0.355	0	    0	    0   	-0.247	0	    0.018	
100	0.601	0	    0	    0	    -0.372	0	    -0.020
110	0.651	1.29	0	    0	    -0.582	0	    -0.031	
120	0.564	2.54	0	    0	    -0.748	0	    -0.024	
130	-0.142	3.58	0	    0.047	-0.700	0	    -0.028	
140	-0.677	3.64	0	    0.069	-0.529	0	    -0.032	
150	-0.723	3.14	0	    0.064	-0.475	0	    -0.032	
160	-2.148	2.56	0	    0.081	0	    1.27	-0.027	
170	-2.707	3.97	-0.175	0.126	0	    1.81	0	
180	-2.529	3.76	-0.174	0.128	0    	1.55	0	      ];

% CY_data = [gamma_r B0	B1	B2	B3	B4	B5	B6]
CY_data = [...
0   0       0       0       0       0       0       0       
10	0.096	0.22	0	    0	    0   	0       0   	
20	0.176	0.71	0	    0	    0	    0	    0   	
30	0.225	1.38	0	    0.023	0	    -0.29	0   	
40	0.329	1.82	0	    0.043   0   	-0.59	0   	
50	1.164	1.26	0.121	0	    -0.242	-0.95	0   	
60	1.163	0.96	0.101	0	    -0.177	-0.88	0   	
70	0.916	0.53	0.069	0	    0   	-0.65	0   	
80	0.844	0.55	0.082	0	    0   	-0.54	0   	
90	0.889	0	    0.138	0	    0   	-0.66	0   
100	0.799	0	    0.155	0	    0    	-0.55	0   	
110	0.797	0	    0.151	0	    0	    -0.55	0   	
120	0.996	0	    0.184	0	    -0.212	-0.66	0.34	
130	1.014	0	    0.191	0	    -0.280	-0.69	0.44	
140	0.784	0	    0.166	0	    -0.209	-0.53	0.38	
150	0.536	0	    0.176	-0.029	-0.163	0	    0.27	
160	0.251	0	    0.106	-0.022	0	    0	    0	    
170	0.125	0	    0.046	-0.012	0	    0	    0	    
180 0       0       0       0       0       0       0     ];

% CN_data = [gamma_r C0	C1	C2	C3	C4	C5]
CN_data = [...
0   0       0       0       0       0       0       
10	0.0596	0.061	0	    0	    0	    -0.074	
20	0.1106	0.204	0	    0	    0	    -0.170	
30	0.2258	0.245	0   	0	    0	    -0.380	
40	0.2017	0.457	0	    0.0067	0	    -0.472	
50	0.1759	0.573	0	    0.0118	0	    -0.523	
60	0.1925	0.480	0	    0.0115	0	    -0.546	
70	0.2133	0.315	0	    0.0081	0	     -0.526	
80	0.1827	0.254	0	    0.0053	0	    -0.443	
90	0.2627	0	    0	    0	    0	    -0.508	
100	0.2102  0	    -0.0195	0	    0.0335	-0.492	
110	0.1567	0	    -0.0258	0	    0.0497	-0.457
120	0.0801	0	    -0.0311	0	    0.0740	-0.396	
130	-0.0189	0	    -0.0488	0.0101	0.1128	-0.420	
140	0.0256	0	    -0.0422	0.0100	0.0889	-0.463
150	0.0552	0	    -0.0381	0.0109	0.0689	-0.476
160	0.0881	0	    -0.0306	0.0091	0.0366	-0.415
170	0.0851	0	    -0.0122	0.0025	0   	-0.220	
180 0       0       0       0       0       0       ];

                        
% interpolate in the tables		
A0 = interp1(CX_data(:,1),CX_data(:,2),gamma_r);
A1 = interp1(CX_data(:,1),CX_data(:,3),gamma_r);
A2 = interp1(CX_data(:,1),CX_data(:,4),gamma_r);
A3 = interp1(CX_data(:,1),CX_data(:,5),gamma_r);
A4 = interp1(CX_data(:,1),CX_data(:,6),gamma_r);
A5 = interp1(CX_data(:,1),CX_data(:,7),gamma_r);
A6 = interp1(CX_data(:,1),CX_data(:,8),gamma_r);

B0 = interp1(CY_data(:,1),CY_data(:,2),gamma_r);
B1 = interp1(CY_data(:,1),CY_data(:,3),gamma_r);
B2 = interp1(CY_data(:,1),CY_data(:,4),gamma_r);
B3 = interp1(CY_data(:,1),CY_data(:,5),gamma_r);
B4 = interp1(CY_data(:,1),CY_data(:,6),gamma_r);
B5 = interp1(CY_data(:,1),CY_data(:,7),gamma_r);
B6 = interp1(CY_data(:,1),CY_data(:,8),gamma_r);

C0 = interp1(CN_data(:,1),CN_data(:,2),gamma_r);
C1 = interp1(CN_data(:,1),CN_data(:,3),gamma_r);
C2 = interp1(CN_data(:,1),CN_data(:,4),gamma_r);
C3 = interp1(CN_data(:,1),CN_data(:,5),gamma_r);
C4 = interp1(CN_data(:,1),CN_data(:,6),gamma_r);
C5 = interp1(CN_data(:,1),CN_data(:,7),gamma_r);

% wind coeffisients
CX =  -(A0 + A1*2*ALw/Loa^2 + A2*2*AFw/B^2 + A3*(Loa/B) + A4*(S/Loa) + A5*(C/Loa) + A6*M);
CY =    B0 + B1*2*ALw/Loa^2 + B2*2*AFw/B^2 + B3*(Loa/B) + B4*(S/Loa) + B5*(C/Loa) + B6*A_SS/ALw;
CN =    C0 + C1*2*ALw/Loa^2 + C2*2*AFw/B^2 + C3*(Loa/B) + C4*(S/Loa) + C5*(C/Loa);

% wind forces and moment
tauX = 0.5*CX*rho_a*V_r^2*AFw;
tauY = 0.5*CY*rho_a*V_r^2*ALw;
tauN = 0.5*CN*rho_a*V_r^2*ALw*Loa;

tau_w = [tauX,tauY,tauN]';