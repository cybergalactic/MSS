function xnext = euler2(xdot,x,h)
% EULER2  Integrate a system of ordinary differential equations using 
%	  Euler's 2nd-order method.
%
% xnext = euler2(xdot,x,h)
%
% xdot  - dx/dt(k) = f(x(k)
% x     - x(k)
% xnext - x(k+1)
% h     - step size
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 

xnext = x + h*xdot;