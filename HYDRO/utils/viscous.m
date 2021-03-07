function Bv = viscous(vessel)
%
% [Bv] = viscous(vessel,surge_estimate) computes the viscous damping
%     matrix Bv. 
% 
%     DOF  1,2,6 : Bvii = beta_i * max(Bii(w))  
%
%     DOF 4      : Bv44  = b * B44(w4) where b is the roll 
%                          amplification factor at w4 satisfying:
%
%                  B44_total(w4) = B44(w4) + B44v(w4) 
%                                = 2 * zeta4 * w4 * (Ix + A44(w4))
% where
%
%     beta_i: percentage of max potential damping (typically 0.1-0.2)
%     zeta4:  relative damping ratio at the roll natural frequency 
%
% Inputs:
%    vessel:       MSS vessel structure
%
% Outputs:
%    Bv(1:6,1:6):  Viscous damping matrix
%         
% Author:    Thor I. Fossen
% Date:      2005-05-26
% Revisions: 2021-03-07 Major improvements, new formulas for viscous damping

[dof1,dof2,nfreqsB,nspeeds] = size(vessel.B);  % potential damping dimensions
w      = vessel.freqs;
nfreqs = length(w);                            % number of freqs for Bv

disp(' ')
disp('-------------------------------------------------------------------')
disp('VISCOUS DAMPING')
disp('-------------------------------------------------------------------')
disp('Viscous damping formula in surge, sway and yaw: Bv(w) = k * max(B(w))')   
visc_input1 = input('Damping factors at w = 0 for DOFs 1,2,6 (default: k = [0.2 0.2 0.2]) : ');

% test for default values
if isempty(visc_input1)
    beta_visc = [0.2 0.2 0.2];
else
    beta_visc = visc_input1;    
end

disp(' ');
disp('Roll damping: B44_total(w4) = B44(w4) + B44v(w4) at w4.')
disp('The relative damping ratio zeta4 is specified such that:')
disp('B44_total(w4) = 2 * zeta4 * w4 * (Ix + A44(w4)).') 

winf = length(w);  % inf frequency

% period DOF 4
MRB = vessel.MRB;
for k = 1:nfreqs
    A(:,:,k) = reshape(vessel.A(:,:,k,1),6,6,1);
end
G   = reshape(vessel.C(:,:,nfreqs,1),6,6);
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

% increase roll damping
disp(' ')
fprintf('Relative damping ratio in roll at w4 = %3.2f rads/s is zeta4 = %5.4f\n',w4,zeta(4));   
fprintf('This can be increased to %3.2f by adding viscous damping B44v = %0.3g\n',zeta_new,B44v_def);
disp(' ')
fprintf('Potential damping: B44(w4) = %0.3g\n',B44);
visc_input2 = input(sprintf('Additional viscous ROLL damping (default: B44v(w4) = %0.3g)       : ',B44v_def));

% test for default values
if isempty(visc_input2)
    B44v = B44v_def;  % default
else
    B44v = visc_input2;
end

% viscous damping factors in surge, sway and yaw
B11_max = beta_visc(1)*max(reshape(vessel.B(1,1,:,1),1,nfreqsB));
B22_max = beta_visc(2)*max(reshape(vessel.B(2,2,:,1),1,nfreqsB));
B66_max = beta_visc(3)*max(reshape(vessel.B(6,6,:,1),1,nfreqsB));

% build diagonal viscous damping matrix for all speeds
Bv =  diag([B11_max B22_max 0 B44v 0 B66_max]);

% compute time constants/periods with and without viscous damping Bv
vessel_new = vessel; 
vessel_new.Bv = Bv;
disp(' ');
disp('VISCOUS DAMPING:')
DPperiods(vessel_new,1);

vessel_new.Bv = 0*Bv;
disp('DAMPING FROM HYDRODYNAMIC CODE:')
DPperiods(vessel_new,1);

disp(' ');
