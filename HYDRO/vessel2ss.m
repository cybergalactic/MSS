function vesselABC = vessel2ss(vessel)
% VESSEL2SS (MSS Hydro)
%
% vesselABC = vessel2ss(vessel,order44,Bvflag) computes the
% hydrodynamic coefficients, retardation function transfer functions using
% SI. The vessel data structure must be generated using ShipX (VERES) or WAMIT. 
%
% Load or generate MSS vessel structure:
%
%   >>load('myship')                     - load vessel struture from myship.mat
%   >>vessel = veres2vessel('s175')      - compute vessel structure    
%   >>vessel = veres2vessel('supply')    - compute vessel structure 
%   >>vessel = wamit2vessel('tanker')    - compute vessel structure 
%
% Generate MSS fluid memory state-space model:
%
%   >>vesselABC = vessel2ss(vessel)
%
% The computed data are stored in: *ABC.mat where * is equal to
% vessel.main.name. The data file *ABC.mat is used by the Simulink templates
%
% Retardation functions and corresponding state-space models are computed
% using potential coefficients with optional viscous damping terms B11v,
% B22v, B44v and B66v. These are computed in viscous_skin.m.
%
% Inputs:
% vessel      MSS vessel structure
%
% Outputs the structure vesselABC with
%   G           = spring stiffness
%   MRB         = rigid-body mass matrix
%   MA          = zero frequency added mass matrix
%   Minv        = inv(MRB+MA)
%   Ar,Br,cr,Dr = zero speed state space models for fluid memory effect
%   r_g         = [xg,yg,zg] vector from O to G        
%
% Author:     Thor I. Fossen
% Date:       2004-08-25  
% Revisions:  2008-02-15  Minor bug fixes
%             2009-09-10  Frequency-domain implementation of fluid
%                         memory effects using the Matlab FDI tools  
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%%
close all

%% ------------------------------------------------------------------------
% USER INPUTS
% -------------------------------------------------------------------------
rho_w = 1025;      % density of water (kg/ms^2)

% system identification options
FDIopt.OrdMax     = 20;
FDIopt.Method     = 2;
FDIopt.Iterations = 20;
FDIopt.LogLin     = 1;
FDIopt.wsFactor   = 0.1;
FDIopt.wminFactor = 0.1;
FDIopt.wmaxFactor = 5;
                                    
%--------------------------------------------------------------------------
%% COMPUTATIONS (do not edit below this line unless you know what you are doing)
% -------------------------------------------------------------------------
disp(' ');
disp('************************* MSS Hydro ********************************')
disp('vesselABC = vessel2ss(vessel) computes the hydrodynamic coefficients')
disp('and retardation functions approximated by transfer functions (fluid')
disp('memory effects) from the MSS <vessel> structure. The results are stored');
disp('in the data file: <vessel name>ABC.mat.')
disp(' ');
disp('Author: Thor I. Fossen');
disp('**********************************************************************')

% check to see if the data are Veres (2D) or WAMIT data (3D)
WAMIT = 0;
if vessel.freqs(end) == 10
    WAMIT = 1;
end

% check the range of vessel.freqs
if vessel.freqs(1)==0 & vessel.freqs(end)~=10
    error('The WAMIT frequencies must be defined such that: vessel.freqs(1:N)= [0 10] where 10 corresponds to the infinite frequency')
end
if vessel.freqs(end)==10 & vessel.freqs(1)~=0
    error('The WAMIT frequencies must be defined such that: vessel.freqs(1:N)= [0 10] where 10 corresponds to the infinite frequency')
end
if vessel.freqs(1)<0.1 & vessel.freqs(end)~=10
    error('The VERES frequencies must be defined such that: vessel.freqs(1:N)= [wmin wmax] where wmin > 0.1 and wmax < 10')
end
if vessel.freqs(1)<0 | vessel.freqs(end)>10
    error('The frequencies must be defined such that: vessel.freqs(1:N)= [wmin wmax] where wmin >= 0 and wmax <= 10')
end

%--------------------------------------------------------------------------
%% Frequency-domain identification
%--------------------------------------------------------------------------
Nmax   = length(vessel.freqs);

if WAMIT == 1
    Nf = Nmax-1;    % for WAMIT computations, remove infinite frequency data
else
    Nf = Nmax;      % for Veres use all frequencies, no infinite frequency data
end

 % frequencies
w = vessel.freqs(1:Nf)';   % does not include ininite frequency

% rigid-body data
MA     = reshape(vessel.A(:,:,Nmax,1),6,6); 
MRB    = vessel.MRB;
G      = reshape(vessel.C(:,:,Nmax,1),6,6);
M      = MRB + MA;

LCG   = vessel.main.CG(1);
VCG   = vessel.main.CG(3);
T_WL  = vessel.main.T;
r_g = [LCG 0 T_WL-VCG];

% indeces for starboard-port symmetric vessels
if WAMIT == 1 % WAMIT
    idx_M = {1,1,1,2,2,2,3,3,4,4,5,6};
    idx_N = {1,3,5,2,4,6,3,5,4,6,5,6};
else % VERES (only element 1,1, in surge, no coupling)
    idx_M = {1,2,2,2,3,3,4,4,5,6};
    idx_N = {1,2,4,6,3,5,4,6,5,6};
end
    
for i = 1:length(idx_M),
    
    dof =[idx_M{i},idx_N{i}];   % 12 elements due to starboard-port symmetry
    
    % hydrodynamic raw data
    A    = reshape( vessel.A(dof(1),dof(2),1:Nf),1,Nf)';
    B    = reshape( vessel.B(dof(1),dof(2),1:Nf),1,Nf)';  
      
    % Call FDIRadMod (Matlab tool for frequency identification)
    if WAMIT == 0   % VERES: identify the infinite-frequency added mass
        
        FDIopt.AinfFlag = 0;
        [Krad,Ainf_hat] = FDIRadMod(w,A,0,B,FDIopt,dof);
        
        MA(dof(1),dof(2)) = Ainf_hat; % estimated infinite added mass matrix
        
    elseif WAMIT == 1 % WAMIT: use the infinite-frequency added mass
        
        Ainf = vessel.A(dof(1),dof(2),Nmax);
        
        if sum(vessel.A(dof(1),dof(2),:)) ~= 0
            FDIopt.AinfFlag = 1;
            [Krad,Ainf_hat] = FDIRadMod(w,A,Ainf,B,FDIopt,dof);
        end
                   
    end
    
    % compute state-space model (Ar,Br,Cr)
    if sum(vessel.A(dof(1),dof(2),:)) == 0
        Arad = 0;
        Brad = 0;
        Crad = 0;
    else
        [Arad,Brad,Crad,Drad] = tf2ss(Krad.num{1},Krad.den{1});
    end
    
    Ar{dof(1),dof(2)} = Arad;
    Br{dof(1),dof(2)} = Brad;
    Cr{dof(1),dof(2)} = Crad;
    Dr{dof(1),dof(2)} = 0;

end

% no data for 1,3 and 1,5 in, only 1,1
if WAMIT == 0   
    Ar{1,3} = 0;
    Br{1,3} = 0;
    Cr{1,3} = 0;
    Dr{1,3} = 0;
    
    Ar{1,5} = 0;
    Br{1,5} = 0;
    Cr{1,5} = 0;
    Dr{1,5} = 0;    
end

% symmetric elements (WAMIT and VERES)
Ar{3,1} = Ar{1,3};
Br{3,1} = Br{1,3};
Cr{3,1} = Cr{1,3};
Dr{3,1} = Dr{1,3};

Ar{4,2} = Ar{2,4};
Br{4,2} = Br{2,4};
Cr{4,2} = Cr{2,4};
Dr{4,2} = Dr{2,4};

Ar{5,1} = Ar{1,5};
Br{5,1} = Br{1,5};
Cr{5,1} = Cr{1,5};
Dr{5,1} = Dr{1,5};

Ar{5,3} = Ar{3,5};
Br{5,3} = Br{3,5};
Cr{5,3} = Cr{3,5};
Dr{5,3} = Dr{3,5};

Ar{6,2} = Ar{2,6};
Br{6,2} = Br{2,6};
Cr{6,2} = Cr{2,6};
Dr{6,2} = Dr{2,6};

Ar{6,4} = Ar{4,6};
Br{6,4} = Br{4,6};
Cr{6,4} = Cr{4,6};
Dr{6,4} = Dr{4,6};

%--------------------------------------------------------------------------
%% Save data
% -------------------------------------------------------------------------
vesselABC.Ar        = Ar;
vesselABC.Br        = Br;
vesselABC.Cr        = Cr;
vesselABC.Dr        = Dr;

vesselABC.MRB       = MRB;
vesselABC.MA        = MA;
vesselABC.G         = G;
vesselABC.Minv      = inv(M);
vesselABC.r_g       = r_g;

save([vessel.main.name 'ABC'],'vesselABC');
disp('------------------------------------------------------------------ ')
disp(['Simulink structure <vesselABC> has been saved to: ' vessel.main.name 'ABC.mat']);
disp('------------------------------------------------------------------ ')
