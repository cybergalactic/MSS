function Bv = viscous(vessel)
%
% [Bv] = viscous(vessel,surge_estimate) computes the viscous damping
%     matrix Bv. 
% 
%     DOF  1,2,6 : Bvii = beta_i * exp(- alpha * w)
%                         where alpha is the exponential rate       
%
%     DOF 4      : Bv44  = b * B44(w4) where b is the roll 
%                          amplification factor at w4 satisfying:
%
%                  B44_total(w4) = B44(w4) + B44v(w4) 
%                                = 2 * zeta4 * w4 * (Ix + A44(w4))
% Inputs:
%    vessel:       MSS vessel structure
%
% Outputs:
%    Bv(:,:,w):    Viscous damping matrix
%         
% Author:    Thor I. Fossen
% Date:      2005-05-26
% Revisions: 2021-03-08 Major improvements, new formulas for viscous damping

w      = vessel.freqs;
Nfreqs = length(w);           % number of freqs for Bv

MRB = vessel.MRB;
for k = 1:Nfreqs
    A(:,:,k) = reshape(vessel.A(:,:,k,1),6,6,1);   
end
G   = reshape(vessel.C(:,:,Nfreqs,1),6,6);      

disp(' ')
disp('-------------------------------------------------------------------')
disp('VISCOUS DAMPING')
disp('-------------------------------------------------------------------')
disp('Viscous damping in surge, sway and yaw: Bvii(w) = beta_i * exp(- alpha * w)')   
T_input = input('Specify the time constants: Ti = (MRBii+Aii(0))/Bvii(0) for DOFs 1,2,6 (default: [50 100 20]) : ');

alpha = 1;

if isempty(T_input)             % test for default values
    T_input = [50 100 20]; 
end

beta_i = [ ( MRB(1,1) + A(1,1,1) ) / T_input(1) 
           ( MRB(2,2) + A(2,2,1) ) / T_input(2)
           ( MRB(6,6) + A(6,6,1) ) / T_input(3)] ;

disp(' ');
disp('Roll damping: B44_total(w4) = B44(w4) + B44v(w4) at w4.')
disp('The relative damping ratio zeta4 is specified such that:')
disp('B44_total(w4) = 2 * zeta4 * w4 * (Ix + A44(w4)).') 

% Natural period for DOF = 4
w_0 = sqrt(G(4,4)/(MRB(4,4)+A(4,4))); 
w4 = natfrequency(vessel,4,w_0,1);
 
A44 = interp1(vessel.freqs,reshape(vessel.A(4,4,:,1),1,length(vessel.freqs)),w4);
B44 = interp1(vessel.freqs,reshape(vessel.B(4,4,:,1),1,length(vessel.freqs)),w4);

% Total damping zeta_new includes viscous damping B44v
M44 = vessel.MRB(4,4) + A44;

% Guideline to increase roll damping
zeta(4) = 0.5 * B44/M44 * (1/w4);

if zeta(4) < 0.05
   zeta_new = 0.05;
   B44v_def = 2*(zeta_new-zeta(4))*w4*M44;
elseif zeta(4) < 0.10
   zeta_new = 0.10;  
   B44v_def = 2*(zeta_new-zeta(4))*w4*M44;  
elseif zeta(4) < 0.15
   zeta_new = 0.15;    
   B44v_def = 2*(zeta_new-zeta(4))*w4*M44;
elseif zeta(4) < 0.20
   zeta_new = 0.20;    
   B44v_def = 2*(zeta_new-zeta(4))*w4*M44;
elseif zeta(4) < 0.25
   zeta_new = 0.25;    
   B44v_def = 2*(zeta_new-zeta(4))*w4*M44;   
else
   zeta_new = zeta(4);    
   B44v_def = 0;  
end

% Increase roll damping
disp(' ')
fprintf('Relative damping ratio in roll at w4 = %3.2f rads/s is zeta4 = %5.4f\n',w4,zeta(4));   
fprintf('This can be increased to %3.2f by adding viscous damping B44v = %0.3g\n',zeta_new,B44v_def);
disp(' ')
fprintf('Potential damping: B44(w4) = %0.3g\n',B44);
visc_input = input(sprintf('Additional viscous ROLL damping (default: B44v(w4) = %0.3g)       : ',B44v_def));

% Test for default values
if isempty(visc_input)
    Bv44 = B44v_def;  % default
else
    Bv44 = visc_input;
end

% Build diagonal viscous damping matrix 
for k = 1:Nfreqs
    Bv(:,:,k) = zeros(6,6);
    Bv(1,1,k) = beta_i(1) * exp( -alpha * w(k) );   
    Bv(2,2,k) = beta_i(2) * exp( -alpha * w(k) );  
    Bv(4,4,k) = Bv44;
    Bv(6,6,k) = beta_i(3) * exp( -alpha * w(k) );       
end

% Compute time constants/periods with and without viscous damping Bv
vessel_new = vessel; 
vessel_new.Bv = Bv;
disp(' ');
disp('VISCOUS DAMPING:')
plotBv(vessel_new);
DPperiods(vessel_new,1);

vessel_new.Bv = 0*Bv;
disp('DAMPING FROM HYDRODYNAMIC CODE:')
DPperiods(vessel_new,1);

disp(' ');
