function [M,D,G,MRB,MA,rCG,rCB,T_z,T_phi,T_theta] = rig
% [M,D,G,MRB,MA,rCG,rCB,T_z,T_phi,T_theta] = rig
% Computes the 6-DOF model parameters of a semi-submersible including the
% natural periods in heave, roll and pitch (T_z, T_phi, T_theta). 
%
%          Inputs: none   
%
%          Outputs: 6x6 model matrices M, D and G
%                    .
%                  M v + D v + G n = tau   
%
%          as well as MRB and MA. Center of gravivty (rCG) and center of 
%          buoynacy (rCB) are given with respect to CO.        
% 
% Author:  2000-03-12 Thor I. Fossen

data_rig

T_draft = abs(draft);                    % draft

N_l = DATCONF.N_l;

rho = DATCONF.rho;
g   = DATCONF.g;
m   = DATCONF.m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for pontoons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main diemensions
H_p  = DATCONF.H_p;            % pontoon height
B_p  = DATCONF.B_p;            % pontoon width
L_p  = DATCONF.L_p;            % pontoon length
CB_p = DATCONF.CB_p;           % pontoon block coefficient

% volume centers in COO
R_p(1,1:3) = [DATCONF.R_p(1,1) DATCONF.R_p(1,2) -H_p/2 ];
R_p(2,1:3) = [DATCONF.R_p(2,1) DATCONF.R_p(2,2) -H_p/2 ];

% deplasement of pontoons
nabla_p = 2*CB_p*(L_p*B_p*H_p);
m_p = rho*nabla_p;

% added mass in surge

alpha_p = 0.10;   % ratio: added mass/mass in surge

% strip theory 
% 2D box shaped data pp. 56-57, Trinatafyllou & Amzallag, MITSG 85-30TN
ratio1   = [0 0.5  1.0  2.0  3.0  4.0  5.0 ];                 % data Figure 2.2
Cm_tab   = [1 1.32 1.51 1.69 1.82 1.91 2.00];

ratio2   = [0 0.125 0.25 0.375 0.50 0.625 0.75 0.875 1.00 ];  % data Figure 2.3
Cr_tab   = [1.0 1.20  1.20 1.20  1.20 1.20  1.26 1.47  1.86 ];

% sway: Cm(Bp/Hp) 
BpdivHp = B_p/H_p;
if BpdivHp>5.0, BpdivHp=5.0; disp('WARNING: Cm22max is used for Cm22'); end
Cm22   = interp1( ratio1, Cm_tab, BpdivHp );
A22_2D = Cm22*rho*pi*(H_p/2)^2;

% heave: Cm(Hp/Bp)
HpdivBp = H_p/B_p;
if HpdivBp>5.0, HpdivBp=5.0; disp('WARNING: Cm33max is used for Cm33'); end
Cm33   = interp1( ratio1, Cm_tab, HpdivBp );
A33_2D = Cm33*rho*pi*(B_p/2)^2;

% roll: Cr(Hp/Bp) 
HpdivBp = H_p/B_p;
if HpdivBp>1.0, HpdivBp=1.0; disp('WARNING: Crmax is used for Cr'); end
Cr =  interp1(ratio2,Cr_tab,HpdivBp);

% show interpolations graphically
subplot(311)
x  = 0:0.1:5;
y  = interp1(ratio1,Cm_tab,x);
plot(ratio1,Cm_tab,'o',x,y,'g',B_p/H_p,Cm22,'sr'),grid,title('C_{m22}')

subplot(312)
x  = 0:0.1:5;
y  = interp1(ratio1,Cm_tab,x);
plot(ratio1,Cm_tab,'o',x,y,'g',H_p/B_p,Cm33,'sr'),grid,title('C_{m33}')

subplot(313)
x  = 0:0.05:1.25;
y  = interp1(ratio2,Cr_tab,x);
plot(ratio2,Cr_tab,'o',x,y,'g',H_p/B_p,Cr,'sr'),grid,title('C_r')

% reduction of added mass in heave due to legs on the top of the pontoons: A33_cyl = 0
k33     = 0.7;    

% added mass with respect to pontoon volume centers
A11 = alpha_p*m_p;
A22 = A22_2D*L_p;
A33 = k33*A33_2D*L_p;
A44 = (1/4)*Cr*pi*rho*L_p*(B_p/2)^2;  
A55 = A33_2D*(2/3)*(L_p/2)^3;
A66 = A22_2D*(2/3)*(L_p/2)^3;

A = diag([A11 A22 A33 A44 A55 A66]);

% transform added mass from pontoon volume centers to COO
A_p1 = Hmtrx(R_p(1,:))'*A*Hmtrx(R_p(1,:));
A_p2 = Hmtrx(R_p(2,:))'*A*Hmtrx(R_p(2,:));

A_p = A_p1 + A_p2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for legs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% volume centers

for i=1:N_l
   R_l(i,:) = [DATCONF.R_l(i,1) DATCONF.R_l(i,2) -T_draft/2 ];   
end

% main dimensions

H_l = T_draft-H_p;                   % heigh of submerged part of legs
CB_l = DATCONF.CB_l;                 % leg block coefficient

for i=1:N_l
   B_l(i) = DATCONF.B_l(i);          % width  of leg 
   L_l(i) = DATCONF.L_l(i);          % length of leg 
   
   if DATCONF.legs == 'box'
      nabla_li(i) = CB_l*B_l(i)*H_l*L_l(i);  % box deplasements
   else % cylinder
      diameter(i) = B_l(i);
      nabla_li(i)  = H_l*pi*(diameter(i)/2)^2;  % cylinder deplasements 
   end
end

% volume centers in COO
for i=1:N_l
   R_l(i,:) = [DATCONF.R_l(i,1) DATCONF.R_l(i,2) -(H_p+H_l/2) ];
end

% added mass for N_l legs
A_l = zeros(6,6);

for i=1:N_l
   if DATCONF.legs == 'box'
      % sway: Cm(Bp/Lp) 
      Cm22   = interp1( ratio1, Cm_tab, B_l(i)/L_l(i) );
      A22_2D = Cm22*rho*pi*(L_l(i)/2)^2;

      % surge: Cm(Lp/Bp)
      Cm11   = interp1( ratio1, Cm_tab, L_l(i)/B_l(i) );
      A11_2D = Cm11*rho*pi*(B_l(i)/2)^2;

      A11 = A11_2D*H_l;
      A22 = A22_2D*H_l;
      A33 = 0;
      A44 = A11_2D*(2/3)*(H_l/2)^3; 
      A55 = A22_2D*(2/3)*(H_l/2)^3;
      A66 = 0;
   else 
      % smooth cylinder Morrisons equation
      Cm = 2.0;                 

      A11 = Cm*rho*(diameter(i)/2)^2;
      A22 = A11;
      A33 = 0;
      A44 = (A11/H_l)*(2/3)*(H_l/2)^3;
      A55 = A44;
      A66 = 0;   
   end

   % added mass for leg i with respect to volume center
   A = diag([A11 A22 A33 A44 A55 A66]);

   % transform added mass from volume center to COO
   A_l = A_l + Hmtrx(R_l(i,:))'*A*Hmtrx(R_l(i,:));
   
end

% total deplacement
nabla_l = sum(nabla_li);
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for bracings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_b = DATCONF.D_b;      % average bracing diameter
L_b = DATCONF.L_b;      % total length of bracing
R_b = DATCONF.R_b;      % bracing volume center coordinates

nabla_b = L_b*pi*(D_b/2)^2;
m_b =rho*nabla_b;

A11 = 1.0*m_b;
A22 = 1.0*m_b;
A33 = 0.1*m_b;

A_b = diag([A11 A22 A33 0 0 0]);

% transform added mass from volume center to COO
A_b = Hmtrx(R_b)'*A_b*Hmtrx(R_b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rigid-body mass/deplasement
nabla  = nabla_p + nabla_l + nabla_b;  % total deplasement
m_disp = rho*nabla;                    % mass of displaced water

disp(sprintf('\n \n \n'))
disp(sprintf('%s %4.1f %s','Draft (COO) =',-T_draft,'m'))
disp(sprintf('%s %4.0f %s','Mass of the vessel =',m,'kg'))
disp(sprintf('%s %4.0f %s \n','Mass of displaced water =',m_disp,'kg'))

Delta_nabla = (m-m_disp)/(rho);
disp(sprintf('%s %4.0f %s','Mass error =',m-m_disp,'kg')) 
disp(sprintf('%s %4.0f %s \n','Deplasement error =',Delta_nabla,'m^3')) 

Delta_Db = 2*sqrt( (nabla_b+Delta_nabla) / (L_b*pi))-D_b;
Delta_Lb = Delta_nabla/( pi*(D_b/2)^2 );

A_legs = sum(B_l.*L_l);
Delta_draft = Delta_nabla/(CB_l*A_legs); 

disp(sprintf('%s %4.3f %s','Alt. 1: change bracing diameter with =',Delta_Db,'m'))
disp(sprintf('%s %4.3f %s','Alt. 2: change length of bracings with =',Delta_Lb,'m'))
disp(sprintf('%s %4.3f %s','Alt. 3: change draft with =',Delta_draft,'m'))
disp(sprintf('%s %4.0f %s \n\n','Alt. 4: change mass with =',m_disp-m,'kg'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rigid Body Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ix    = m*DATCONF.r_x^2;                 % moments of inertia
Iy    = m*DATCONF.r_y^2;                 
Iz    = m*DATCONF.r_z^2;                 

% water plane area

if DATCONF.legs == 'box'
   Ao = CB_l*(B_l.*L_l);       % areas of box shaped legs
else
   Ao = pi*(0.5*B_l).^2;       % areas of cylinder shaped legs
end

Awp   = sum(Ao)               % water plane area

% Moment of areas
for i=1:N_l
   if DATCONF.legs == 'box'
      I_Li(i) = Ao(i)*( R_l(i,1)^2 + (1/12)*L_l(i)*B_l(i)^3 );    
      I_Ti(i) = Ao(i)*( R_l(i,2)^2 + (1/12)*B_l(i)*L_l(i)^3 );    
   else % cylinders
      I_Li(i) = Ao(i)*( R_l(i,1)^2 + (1/4)*Ao(i)/pi ); 
      I_Ti(i) = Ao(i)*( R_l(i,2)^2 + (1/4)*Ao(i)/pi );          
   end
end

I_L = sum(I_Li);
I_T = sum(I_Ti);

% Center of buoynacy and WP area center in COO
sum_l = zeros(1,3);
sum_wp = zeros(1,3);
for i=1:N_l
   sum_l  = sum_l + R_l(i,:)*(Ao(i)*H_l);
   sum_wp = sum_wp + [R_l(i,1:2) -T_draft]*Ao(i);
end
sum_p = 0.5*nabla_p*(R_p(1,:)+R_p(2,:));

rCB = (sum_p + sum_l + nabla_b*R_b )/nabla   % CB
rWP = sum_wp/Awp;                            % water plane area center

% Metacentic heights
rCG = DATCONF.R_CG'

BG = rCB(3)-rCG(3);

GM_L = I_L/nabla-BG
GM_T = I_T/nabla-BG

% Stiffness matrix

G33 = rho*g*Awp;
G44 = rho*g*nabla*GM_T;
G55 = rho*g*nabla*GM_L;

% G-matrix with respect to water plane area
G_WP = diag([0 0 G33 G44 G55 0]);

% transform G_Wp to COO
G = Hmtrx(rWP)'*G_WP*Hmtrx(rWP);

% transform G_COO to CG
G_CG = Hmtrx(-rCG)'*G*Hmtrx(-rCG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of M_COS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MA = A_p + A_l + A_b;

% transform from COO to CG
MA_CG  = Hmtrx(-rCG)'*MA*Hmtrx(-rCG);

% rigid body mass
MRB_CG = [ m*eye(3)     zeros(3,3)
           zeros(3,3)  diag([Ix Iy Iz]) ];

MRB = Hmtrx(rCG)'*MRB_CG*Hmtrx(rCG);

% Inertia matrix in COO, COS and CG
M = MA + MRB;
M_CG = MA_CG + MRB_CG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Damping D_COS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_z     = 0.995;    % w_damped = beta * w_undamped 
beta_phi   = 0.995;    % where beta = reduction factor
beta_theta = 0.995;

v_z     = 2*sqrt(1-beta_z^2);
v_phi   = 2*sqrt(1-beta_phi^2);
v_theta = 2*sqrt(1-beta_theta^2);

% damping (diagonal structure at CG)
D_CG11 = M_CG(1,1)/DATCONF.T_x;
D_CG22 = M_CG(2,2)/DATCONF.T_y;
D_CG66 = M_CG(6,6)/DATCONF.T_n;

D_CG33 = v_z    *sqrt(M_CG(3,3)*G_CG(3,3));
D_CG44 = v_phi  *sqrt(M_CG(4,4)*G_CG(4,4));
D_CG55 = v_theta*sqrt(M_CG(5,5)*G_CG(5,5));

D_CG = diag([D_CG11 D_CG22 D_CG33 D_CG44 D_CG55 D_CG66]);

% transform D_CG to COO
D = Hmtrx(rCG)'*D_CG*Hmtrx(rCG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute some useful parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_syst = [ zeros(6,6)      eye(6)
          -inv(M)*G   -inv(M)*D ];
       
eig_syst = eig(A_syst);

figure(2),figure(gcf)
subplot(111)
plot(eig_syst(4:12),'r^'),title('6 DOF eigenvalues'),grid

T_z     = 2*pi/sqrt( (G(3,3)/M(3,3))*(1-(v_z/2)^2)     );
T_phi   = 2*pi/sqrt( (G(4,4)/M(4,4))*(1-(v_phi/2)^2)   );
T_theta = 2*pi/sqrt( (G(5,5)/M(5,5))*(1-(v_theta/2)^2) );

save rig MA MRB D G


