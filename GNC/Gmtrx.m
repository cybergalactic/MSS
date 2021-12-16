function G = Gmtrx(nabla,A_wp,GMT,GML,LCF,r_bp)
% G = Gmtrx(nabla,A_wp,GMT,GML,LCF,r_bp) computes the 6x6 system spring 
% stiffness matrix G about an arbitrarily point P for a floating vessel 
% (small roll and pitch angles). For submerged vessels, use gvect.m or
% gRvect.m.
% 
% Inputs:  
%  nabla: volume displacement (m3)
%  Awp: waterplane area (m2)
%  GMT, GML: transverse/longitudinal metacentric heights (m)
%  LCF = x coordinate from the CO to the CF, positive forwards negative for
%        conventional ships (m)
%  r_bp = [x_p y_p z_p]': location of the point P with respect to the CO (m),  
%          use r_bp = [0, 0, 0]' for CO midships
%
% Author:     Thor I. Fossen
% Date:       14 Jun 2001
% Revisions:  26 Jun 2002 variable Awp was replaced by A_wp, one zero
%                         in G was removed (bugfix)
%             20 Oct 2008 updated the documentation, start using r_p for an
%                         arbitrarily point
%             25 Apr 2019 added LCF as input parameter
%             16 Dec 2021 minor updates of the documentation

rho = 1025;  % density of water
g   = 9.81;	 % acceleration of gravity

% Location of the center of flotation (CF)
r_bf = [LCF, 0, 0]';

% Hydrostatic quantities expressed in the CF
G33_CF  = rho * g * A_wp;
G44_CF  = rho * g * nabla * GMT;
G55_CF  = rho * g * nabla * GML;
G_CF = diag([0 0 G33_CF G44_CF G55_CF 0]);  
G_CO = Hmtrx(r_bf)' * G_CF * Hmtrx(r_bf);
G = Hmtrx(r_bp)' * G_CO * Hmtrx(r_bp);

end