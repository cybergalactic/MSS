function [Anew,Bnew,Bv] = ABCtransform(vessel,striptheory,plot_flag)
% ABCTRANSFORM (MSS Hydro)
%
% [Bnew,Anew,Bv,B_roll]] = ABCtransform(vessel,striptheory,plot_flag)
%
% inputs:  vessel                               MSS vessel structure
%          striptheory = {'veres'}              software program
%          plot_flag                            1 = plots raw data in MSS axes                                                 
%                                               0 or empty = no plotting
% 
% outputs: Anew                                 new A terms
%          Bnew                                 new B terms
%          Bv                                   viscous damping
%          Bv44                                 viscous roll damping
% 
% Author:    Thor I. Fossen
% Date:      2007-08-24
% Revisions: 2009-09-10 new estimate for A11
%                       using new viscous damping function viscous.m
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%%

Nfreqs  = length(vessel.freqs);
Nspeeds = length(vessel.velocities);

% zero speed terms
speed = 1;
A_0    = reshape(vessel.A(:,:,:,speed),6,6,Nfreqs);
B_0    = reshape(vessel.B(:,:,:,speed),6,6,Nfreqs);

% make zero speed matrices symmetric to account for numerical errors
for i=1:Nfreqs,
   A_0(1:6,1:6,i) = 0.5*(A_0(1:6,1:6,i) + A_0(1:6,1:6,i)');
   B_0(1:6,1:6,i) = 0.5*(B_0(1:6,1:6,i) + B_0(1:6,1:6,i)');   
end

% initial value
Bnew = vessel.B;        
Anew = vessel.A;  

if strcmp(striptheory,'veres')

    % ---------------------------------------------------------------------
    % compute: A11(0) = -2.7*(rho/L^2)*nabla^(5/3) c.f. S�ding (1982)
    % ---------------------------------------------------------------------
    nabla = vessel.main.nabla;
    Lpp   = vessel.main.Lpp;
    rho   = vessel.main.rho;
    A11_0 = 2.7*(rho/Lpp^2)*nabla^(5/3);  
    
    A22_0 = vessel.A(2,2,1,1);
    
    % ratio
    alpha = A11_0/A22_0;
    
    % ---------------------------------------------------------------------
    % compute A11(w) and B11(w) by scaling the A22(w) and B22(w)data
    % ---------------------------------------------------------------------
    for speed=1:Nspeeds,
        for i=1:Nfreqs,
            Anew(1,1,i,speed) = alpha*Anew(2,2,i,speed);
            Bnew(1,1,i,speed) = alpha*Bnew(2,2,i,speed);
        end
    end
    
end

%--------------------------------------------------------------------------
%% viscous damping
%--------------------------------------------------------------------------
vessel_new = vessel;
vessel_new.A = Anew;
vessel_new.B = Bnew;

Bv = viscous(vessel_new);

%--------------------------------------------------------------------------
%% plot
%--------------------------------------------------------------------------
if plot_flag == 1,
  plot_speedterms(vessel,Anew,Bnew,Bv)
end
