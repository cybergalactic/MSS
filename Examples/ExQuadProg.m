% ExQuadProg Quadratic programming for way-point trajectory generation.
%            Two cubic polynominals are fitted to two way-points where the 
%            speed and position are specified. The starting time is t0 while
%            the arrival time t1 at next way-point is unknown.
%
%            Cubic spline between two points:
%            x(t) = a3*t^3 + a2*t^2 + a1*t + a0
%            y(t) = b3*t^3 + b2*t^2 + b1*t + b0
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2004 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see http://www.gnu.org/licenses
% 
% E-mail: contact@marinecontrol.org
% URL:    http://www.marinecontrol.org

% Author:    Thor I. Fossen
% Date:      9 July 2002
% Revisions: 

Umax = 10;    % maximum speed

% way-point 0
t0 = 0;       % time
x0 = 10;      % x-position 
y0 = 10;      % y-psoition
U0 = 0;       % speed

% way-point 1 
% time t1 is unknown (will be optimized)
x1 = 200;     % x-position 
y1 = 100;     % y-psoition
U1 = 5;       % speed

% direction from way-point 0 to way-point 1
psi0 = atan2(y1-y0,x1-x0);
psi1 = psi0;                  % last way-point equals direction of previous

% --------------------------------------------------------------------
% MAIN LOOP for computation of the 9 unknowns,
% time t1 and the polynominal parameters x = [a3 a2 a1 a0 b3 b2 b1 b0]'
% --------------------------------------------------------------------

Jbar_min = 1e10;

for dt = 1:100,
    
    t1 = t0+dt;  
        
    y = [x0 x1 y0 y1 U0*cos(psi0) U0*sin(psi0) U1*cos(psi1) U1*sin(psi1)]';
    
    A = [t0^3    t0^2 t0 1 0 0 0 0
        t1^3    t1^2 t1 1 0 0 0 0
        0 0 0 0 t0^3 t0^2 t0 1
        0 0 0 0 t1^3 t1^2 t1 1  
        3*t0^2 2*t0   1  0 0 0 0 0 
        0 0 0 0 3*t0^2 2*t0   1  0      
        3*t1^2 2*t1   1  0 0 0 0 0 
        0 0 0 0 3*t1^2 2*t1   1  0   ];

    H  = A'*A;
    f  = -y'*A;
    
    [x,J] = quadprog(H,f,[],[]);
    
    Jbar = 0.5*(J - y'*y); 
    
    % speed profile between way-points 0 and 1
    t = linspace(t0,t1,101);
    Acoeff = x(1:4);
    Bcoeff = x(5:8);
    xdot = polyval([0 3*Acoeff(1) 2*Acoeff(2) Acoeff(3)],t);
    ydot = polyval([0 3*Bcoeff(1) 2*Bcoeff(2) Bcoeff(3)],t);
    U    = sqrt(xdot.^2+ydot.^2);
    
    % select optimal solution 
    if max(U) < Umax,
        xopt = x;
        topt = t1;
        break
    end
    
end

% optimal solution
t1 = topt
Acoeff = xopt(1:4)  
Bcoeff = xopt(5:8)

% --------------------------------------------------------------------
% graphics
% --------------------------------------------------------------------

t = linspace(t0,t1,101);

figure(gcf)
subplot(411)
plot(t,polyval(Acoeff,t),'linewidth',2);
grid

subplot(412)
plot(t,polyval(Bcoeff,t),'linewidth',2);
grid

subplot(413)
plot([y0 y1],[x0 x1],'ro',polyval(Bcoeff,t),polyval(Acoeff,t),'linewidth',2);
axis('equal'),grid    

subplot(414)
xdot = polyval([0 3*Acoeff(1) 2*Acoeff(2) Acoeff(3)],t);
ydot = polyval([0 3*Bcoeff(1) 2*Bcoeff(2) Bcoeff(3)],t);
U    = sqrt(xdot.^2+ydot.^2);
plot(t,U,'linewidth',2);
grid    