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
% nonlinear model to describe stall:
%
%  CL(alpha) = (1-sigma) * CL_alpha * alpha + ...
%              sigma     * 2 * sign(alpha) * sin(alpha)^2 * cos(alpha) 
%
% where 0 <= sigma <= 1 is a blending parameter. 
%
% Outputs:
%  CL:      Lift coefficient as a function of alpha   
%  CD:      Drag coefficient as a function of alpha   
%
% Inputs:
%  b:       Wing span (m)
%  S:       Wing area (m^2)
%  CD_0:    Parasitic drag (alpha = 0), typically 0.1-0.2 for a streamlined body
%  alpha:   Angle of attack, scalar or vector (rad)
%  sigma:   Blending parameter between 0 and 1, use sigma = 0 for linear lift 
%  display: Use 1 to plot CD and CL (optionally)
%
% Examples:
%
% Cylinder-shaped AUV with length L = 1.8, diameter D = 0.2 and CD_0 = 0.3:
%    b = L; 
%    S = b * D;
%    CD_0 = 0.1;
%    alpha = (-5:1:80)*pi/180;
%    [CL,CD] = coeffLiftDrag(b, S, CD_0, alpha, 0.2, 1)
% 
% Author:    Thor I. Fossen
% Date:      2021-05-25 
%   2024-10-31 Bug fixes and added plot for optimal alpha corresponding to max(CL/CD)

% Oswald efficiency number
% 0.3-0.4 is typically values for a cylinder-shaped AUVs and torpedos
% 0.5-0.7 is typically for streamlined AUVs such as the NPS AUV
e = 0.3; 

% Wing aspect ratio
AR = b^2/S; 

% Linear lift
CL_alpha = pi * AR / ( 1 + sqrt(1 + (AR/2)^2) );
CL_linear = CL_alpha * alpha;  

% Nonlinear lift (blending function)
CL = (1 - sigma) .* CL_linear + ...
    sigma .* 2 .* sign(alpha) .* sin(alpha).^2 .* cos(alpha);

% Parasitic and induced drag
CD = CD_0 + CL.^2 / (pi * e * AR); 

% Optionally plot
if (nargin == 6 && display == 1)
    figure(gcf)
    % Plot Drag Coefficient CD(alpha)
    subplot(311), plot(rad2deg(alpha), CD, 'linewidth', 2)
    title('Drag coefficient CD(\alpha)'), xlabel('\alpha (deg)'), grid
    
    % Plot Lift Coefficient CL(alpha)
    subplot(312), plot(rad2deg(alpha), CL, 'linewidth', 2)
    title('Lift coefficient CL(\alpha)'), xlabel('\alpha (deg)'), grid
    
    % Calculate and Plot Lift-to-Drag Ratio CL/CD
    CL_CD = CL ./ CD;
    
    % Find the index of the maximum CL/CD
    [max_CL_CD, idx_opt] = max(CL_CD);
    alpha_opt = rad2deg(alpha(idx_opt)); 
    
    % Plot CL/CD with optimal alpha indication
    subplot(313)
    plot(rad2deg(alpha), CL_CD, 'linewidth', 2)
    hold on
    plot(alpha_opt, max_CL_CD, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r') 
    title('Lift-to-Drag Ratio CL/CD(\alpha)')
    xlabel('\alpha (deg)'), ylabel('CL/CD')
    legend('CL/CD', ['Optimal \alpha = ', num2str(alpha_opt, '%.2f'), ' deg'])
    grid on
    hold off
end