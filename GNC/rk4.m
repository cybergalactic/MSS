function xnext = rk4(f,x,u,h,t)
% RK4	Integrate a system of ordinary differential equations using
%	    Runge-Kutta's 4th-order method. The control input u can be constant
%	    over the sampling interval h or a time-varying function u = g(x,t).
%
% The output xnext = x(k+1) is:
%
%    xnext = rk4('f',x,'g',h,t)     nonautonomous systems 
%    xnext = rk4('f',x,u,h)         autonomous systems (with constant u)
%
% where
%
%    f     - external function: 
%            dx/dt(k) = f(x(k),u(k),t(k))  nonautonomous systems
%            dx/dt(k) = f(x(k),u(k))       autonomous systems
%
%    x     - x(k)
%    u     - external function: u(k) = g(x(k),t(k))   nonautonomous systems
%            u(k) = constant over the sample time h   autonomous systems
%    h     - step size
%    t     - time t(k) - ONLY NEEDED FOR NONAUTONOMOUS SYSTEMS
%
% Ex 1:   function u = g(x,t)
%           u = -2*x + cos(t);
%
%         function dx = f(x,u,t)
%           dx = -x^2 + u + cos(t);     
%
%         ===>  xnext = rk4('f',x,'g',h,t) 
%
% Ex 2:   u = constant
%
%         function dx = f(x,u), 
%           dx = -x^2 + u;
%
%         ===>  xnext = rk4('f',x,u,h)
%
% Author:    Thor I. Fossen
% Date:      2001-07-14
% Revisions: 2007-11-23 Christian Holden - update for time-varying systems

xo = x;

if nargin == 4   % autonomous sytem

    k1 = h*feval(f,xo,u);
    x  = xo+0.5*k1;
    k2 = h*feval(f,x,u);
    x  = xo+0.5*k2;
    k3 = h*feval(f,x,u);
    x  = xo+k3;
    k4 = h*feval(f,x,u);

else  % nonautonomous sytem with additional input time t
    
    k1 = h*feval(f,xo,feval(u,xo,t),t);
    x  = xo+0.5*k1;
    k2 = h*feval(f,x,feval(u,x,t+h/2),t+h/2);
    x  = xo+0.5*k2;
    k3 = h*feval(f,x,feval(u,x,t+h/2),t+h/2);
    x  = xo+k3;
    k4 = h*feval(f,x,feval(u,x,t+h),t+h);

end

xnext = xo + (k1+2*(k2+k3)+k4)/6;