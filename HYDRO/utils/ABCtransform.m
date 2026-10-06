function [Anew,Bnew] = ABCtransform(vessel,striptheory,plot_flag)
% [Anew,Bnew] = ABCtransform(vessel,striptheory,plot_flag)
%
% inputs:  vessel                               MSS vessel structure
%          striptheory = {'veres'}              software program
%          plot_flag                            1 = plots raw data in MSS axes                                                 
%                                               0 or empty = no plotting
% 
% outputs: Anew                                 new A terms
%          Bnew                                 new B terms
% 
% Author:    Thor I. Fossen
% Date:      2007-08-24
% Revisions: 2009-09-10 new estimate for A11
%            2021-03-07 new estimates for A11 and B11
%            2026-09-28 removed legacy viscous damping computation

Nfreqs  = length(vessel.freqs);
Nspeeds = length(vessel.velocities);

% zero speed terms
speed = 1;
A_0    = reshape(vessel.A(:,:,:,speed),6,6,Nfreqs);
B_0    = reshape(vessel.B(:,:,:,speed),6,6,Nfreqs);

% make zero speed matrices symmetric to account for numerical errors
for i=1:Nfreqs
   A_0(1:6,1:6,i) = 0.5*(A_0(1:6,1:6,i) + A_0(1:6,1:6,i)');
   B_0(1:6,1:6,i) = 0.5*(B_0(1:6,1:6,i) + B_0(1:6,1:6,i)');   
end

% initial value
Bnew = vessel.B;        
Anew = vessel.A;  

if strcmp(striptheory,'veres')

    % ---------------------------------------------------------------------
    % compute: A11(0) = -2.7*(rho/L^2)*nabla^(5/3) c.f. Søding (1982)
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
    for speed=1:Nspeeds
        for i=1:Nfreqs
            Anew(1,1,i,speed) = alpha*Anew(2,2,i,speed);
            Bnew(1,1,i,speed) = alpha*Bnew(2,2,i,speed);
        end
    end
    
end

%--------------------------------------------------------------------------
%% plot
%--------------------------------------------------------------------------
if plot_flag == 1
  plot_speedterms(vessel,Anew,Bnew)
end
