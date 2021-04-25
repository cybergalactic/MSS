function D = Dmtrx(T_126,zeta_45,MRB,MA,hydrostatics)
% D = Dmtrx([T1, T2, T6],[zeta4,zeta5],MRB,MA,hydrostatics)
% computes the 6x6 linear damping matrix for marine craft (submerged and
% floating) by specifying the time constants [T1, T2, T6] in DOFs 1,2 and 6. 
% The time constants can be found by open-loop step responses. For roll and
% pitch the relative damping ratios are specified using [zeta4, zeta5]. 
% For floating vessels it is assumed that zeta3 = 0.2 in heave, while
% submerged vehicles are assumed to be neutrally buoyant, W = B, with equal 
% time constants in heave and sway, that is T3 = T2.
%
% Inputs: T_126 = [T1, T2, T6]: time constants for DOFs 1, 2 and 6
%         zeta_45 = [zeta4, zeta5]: relative damping ratios in DOFs 4 and 5
%         MRB: 6x6 rigid-body system matrix (see rbody.m)
%         MA:  6x6 hydrodynamic added mass system matrix
%         hydrostatics = G for surface craft (see Gmtrx.m)
%         hydrostatics = [W r_bg' r_bb'] for neutrally buoyant submerged 
%           vehicles where W = m*g, r_bg = [xg,yg,zg]' and r_bb = [xb,yb,zb]'
%
% Output: D: 6x6 diagonal linear damping matrix
%
% Examples:  
%  supply: D = Dmtrx([10 50 20], [0.2, 0.3], MRB, MA, G)
%  tanker: D = Dmtrx([20 100 50], [0.2, 0.3], MRB, MA, G)
%  AUV:    D = Dmtrx([5 10 10], [0.2, 0.3], MRB, MA,[m*g, [0,0,zg], [0,0,0])
%
% Author:     Thor I. Fossen 
% Date:       24 Apr 2021
% Revisions:  

M = MRB + MA;

T1 = T_126(1);
T2 = T_126(2);
T6 = T_126(3);
zeta4 = zeta_45(1);
zeta5 = zeta_45(2);

if isvector(hydrostatics) % submerged vehicle: hydrostatics = [W r_bg r_bb]' 
    
    W = hydrostatics(1);       
    r_bg = hydrostatics(2:4);
    r_bb = hydrostatics(5:7);
    
    T3 = T2;       % for AUVs, assume same time constant in sway and heave
    w4 = sqrt( W * (r_bg(3)-r_bb(3)) / M(4,4) );
    w5 = sqrt( W * (r_bg(3)-r_bb(3)) / M(5,5) );

    D = diag( [M(1,1)/T1 M(2,2)/T2 M(3,3)/T3...
        M(4,4)*2*zeta4*w4  M(5,5)*2*zeta5*w5 M(6,6)/T6 ] );   
    
else                      % surface craft: hydrostatics = G
    
    G33 = hydrostatics(3,3);
    G44 = hydrostatics(4,4);
    G55 = hydrostatics(5,5);
    
    zeta3 = 0.2;           % for ships 
    w3 = sqrt( G33 / M(3,3) );
    w4 = sqrt( G44 / M(4,4) );
    w5 = sqrt( G55 / M(5,5) );
    
    D = diag( [M(1,1)/T1 M(2,2)/T2 M(3,3)*2*zeta3*w3...
        M(4,4)*2*zeta4*w4  M(5,5)*2*zeta5*w5 M(6,6)/T6 ] );
    
end
   

       
