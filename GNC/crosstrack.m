function [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y, flag)
%  [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y, flag) computes 
% the coordinate origin (x_p, y_p) of the path-tangential reference frame
% from the reference point (x_ref, y_ref) to the target position (x_p, y_p), 
% and cross-track error y_e for a craft at position (x, y). 
%
% Input:    (x_t,y_t), North-East target positions
%           (x_ref,y_ref). North-East reference positions
%           (x,y), craft North-East positions
%           flag (optionally) = 1 (plot geometry)
%
% Outputs:  (x_p,y_p) origin of path-tangential frame expressed in NED
%           y_e cross-track error
%
% See also crosstrackWpt(x2, y2, x1, y1, x, y), which computes the
% cross-track error between two waypoints (x1,y1) and (x2,y2).
%  
% Author:    Thor I. Fossen
% Date:      10 July 2020
% Revisions: 6 Oct. 2020, Thor I. Fossen, fixed minor typos

% path-tangential angle
pi_p = atan2(y_t-y_ref, x_t-x_ref);

cp = cos(pi_p);
sp = sin(pi_p); 
tp = tan(pi_p);

% A * z = b  -->  x = inv(A) * b
A = [cp sp  0
    -sp cp 1
     tp -1 0 ];

b = [ cp*x + sp*y
     -sp*x + cp*y
      tp*x_t - y_t ];
      
z = inv(A) * b;
x_p = z(1); 
y_p = z(2);
y_e = z(3);

if (nargin == 7 && flag == 1)
    axis normal;
    figure(gcf)
    plot(y,x,'ro',y_t,x_t,'ko',y_p,x_p,'bo',y_ref,x_ref,'ko')
    hold on
    quiver(y_ref,x_ref,y_t-y_ref,x_t-x_ref,0,'-b')
    quiver(y_p,x_p,y-y_p,x-x_p,0,'-r')
    hold off
    grid;
    axis equal;
end
