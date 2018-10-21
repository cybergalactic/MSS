function Bv = viscous(vessel)
% VISCOUS_SKIN (MSS Hydro)
%
% [Bv] = viscous(vessel,surge_estimate) computes a constant viscous skin
% friction and roll damping matirx Bv. 
% 
%     DOF (i = 1,2,6):  Bv_ii = beta_i*max(B_ii(w))  
%
%     DOF 4          :  Bv_44 = b*B44(w4) where b is the roll 
%                               amplification factor at w4 satisfying:
%
%               Bv_44(w4) = (1+b)*B44(w4) = 2*zeta*sqrt(C44*(Ix+A44(w4)))') 
% where:
%
%     beta_i: percentage of max potential damping         (typically 0.1-0.2)
%     zeta:   relative damping ratio at natural frequency (typically 0.05-0.10)
%
% Inputs:
%    vessel.MRB(1:6,1:6):    MSS vessel 
%
% Outputs:
%    Bv(1:6,1:6):  Viscous damping matrix
%         
% Author:    Thor I. Fossen
% Date:      2005-05-26
% Revisions: 
%
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

[dof1,dof2,nfreqsB,nspeeds] = size(vessel.B);  % potential damping dimensions
w      = vessel.freqs;
nfreqs = length(w);                     % number of freqs for Bv

%--------------------------------------------------------------------------
disp(' ')
disp('-----------------------------------------')
disp('VISCOUS DAMPING')
disp('-----------------------------------------')
disp('Skin friction in surge, sway and yaw:')
disp('   B_v(w) = k*max(B_p)  (deault k = 0.2)')   
disp(' ');
disp('Viscous roll damping B44(w4) = B44p(w4) + B44v at resonance frequency w4.')
disp('   The relative damping ratio zeta4 is specified such that:')
disp('   B44(w4) = 2*zeta4*w4*(Ix+A44(w4))') 
winf = length(w);  % inf frequency

% period DOF 4
[T,zeta] = DPperiods(vessel,1);
w4 = 2*pi/T(4);
 
A44 = interp1(vessel.freqs,reshape(vessel.A(4,4,:,1),1,length(vessel.freqs)),w4);
B44 = interp1(vessel.freqs,reshape(vessel.B(4,4,:,1),1,length(vessel.freqs)),w4);

% Total damping zeta_new includes viscous dampng B44v
M44 = vessel.MRB(4,4) + A44;

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
disp(sprintf('Relative damping ratio in roll at resonance w4 = %3.2f rads/s is zeta4 = %5.4f',w4,zeta(4)));   
disp(sprintf('This is increased to %3.2f by adding viscous damping B44v = %0.3g',zeta_new,B44v_def));
disp(' ')
disp(sprintf('Potential damping:        B44p(w4) = %0.3g',B44));
visc_input1 = input(sprintf('Additional viscous ROLL damping (default: B44v(w4) = %0.3g)      : ',B44v_def));
visc_input2 = input('Viscous damping at w = 0 for DOF 1,2,6  (default: k = [0.2 0.2 0.2]) : ');

% test for default values
if isempty(visc_input1)
    B44v = B44v_def;  % default
else
    B44v = visc_input1;
end

if isempty(visc_input2)
    beta_visc = [0.2 0.2 0.2];
else
    beta_visc = visc_input2;    
end

% skin friction factors in surge, sway and yaw
B11_max = beta_visc(1)*max(reshape(vessel.B(1,1,:,1),1,nfreqsB));
B22_max = beta_visc(2)*max(reshape(vessel.B(2,2,:,1),1,nfreqsB));
B66_max = beta_visc(3)*max(reshape(vessel.B(6,6,:,1),1,nfreqsB));

% build diagonal viscous damping matrix for all speeds
Bv =  diag([B11_max B22_max 0 B44v 0 B66_max]);

% compute time constants/periods with viscous damping Bv
vessel_new = vessel;
vessel_new.Bv = Bv;
DPperiods(vessel_new,1);

disp(' ');
