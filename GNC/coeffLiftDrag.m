function [CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,sigma,display)
% [CL,CD] = coeffLiftDrag(b,S,CD_0,alpha,sigma) computes the hydrodynamic 
% lift CL(alpha) and drag CD(alpha) coefficients as a function of alpha
% (angle of attack) of a submerged "wing profile" using (Beard and McLain 2012)
%
%  CD(alpha) = CD_p + (CL_0 + CL_alpha * alpha)^2 / (pi * e * AR)
%  CL(alpha) = CL_0 + CL_alpha * alpha
% 
% where CD_p is the parasitic drag (profile drag of wing, friction and
% pressure drag of control surfaces, hull, etc.), CL_0 is the zero angle 
% of attack lift coefficient, AR = b^2/S is the aspect ratio and e is the  
% Oswald efficiency number. For lift it is assumed that
%
%  CL_0 = 0
%  CL_alpha = pi * AR / ( 1 + sqrt(1 + (AR/2)^2) );
%
% implying that for alpha = 0, CD(0) = CD_0 = CD_p and CL(0) = 0. For
% high angles of attack the linear lift model can be blended with a
% nonlinear model to describe stall
%
%  CL(alpha) = (1-sigma) * CL_alpha * alpha + ...
%              sigma     * 2 * sign(alpha) * sin(alpha)^2 * cos(alpha) 
%
% where 0 <= sigma <= 1 is a blending parameter. 
%
% Outputs:
%  CL:      lift coefficient as a function of alpha   
%  CD:      drag coefficient as a function of alpha   
%
% Inputs:
%  b:       wing span (m)
%  S:       wing area (m^2)
%  CD_0:    parasitic drag (alpha = 0), typically 0.1-0.2 for a streamlined body
%  alpha:   angle of attack, scalar or vector (rad)
%  sigma:   blending parameter between 0 and 1, use sigma = 0 for linear lift 
%  display: use 1 to plot CD and CL (optionally)
%
% Examples:
%
% Cylinder-shaped AUV with length L = 1.8, diameter D = 0.2 and CD_0 = 0.3:
%    alpha = 0.1 * pi/180;
%    [CL,CD] = coeffLiftDrag(0.2, 1.8*0.2, 0.3, alpha, 0.2)
%    alpha = (-5:1:80)*pi/180;
%    [CL,CD] = coeffLiftDrag(0.2, 1.8*0.2, 0.3, alpha, 0.2, 1)
% 
% Author:    Thor I. Fossen
% Date:      25 April 2021 

e = 0.7;             % Oswald efficiency number
AR = b^2/S;          % wing aspect ratio

% linear lift
CL_alpha = pi * AR / ( 1 + sqrt(1 + (AR/2)^2) );
CL = CL_alpha * alpha;  

% parasitic and induced drag
CD = CD_0 + CL.^2 / (pi * e * AR); 

% nonlinear lift (blending function)
CL = (1-sigma) .* CL + sigma .* 2 .* sign(alpha).*sin(alpha).^2.*cos(alpha);

% optionally plot
if (nargin == 6 && display == 1)
    figure(gcf)
    subplot(211),plot(alpha*180/pi,CD,'linewidth',2)
    title('Drag coefficient CD(\alpha)'),xlabel('\alpha (deg)'),grid
    subplot(212),plot(alpha*180/pi,CL,'linewidth',2)
    title('Lift coefficient CL(\alpha)'),xlabel('\alpha (deg)'),grid
end