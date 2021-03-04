% ExPathGen  Path generation using cubic polynominals
% Author:    Thor I. Fossen
% Date:       1 November 2002
% Revisions: 13 February 2012   Replaced Ylabel with ylabel

% Cubic spline between two points
% x(th)     =   a3*th^3 + a2*th^2 + a1*th + a0
% dx/dt(th) = 3*a3*th^2 + 2*a2*th + a1

clear all
% way-point database
wpt.pos.x   = [0 200 400 700 1000];
wpt.pos.y   = [0 200 500 400 1200];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% COMPUTATION OF THE PATH ( x(th) ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select five active waypoints for curve fitting
x1     = wpt.pos.x(1);
x2     = wpt.pos.x(2);    
x3     = wpt.pos.x(3);   
x4     = wpt.pos.x(4);   
x5     = wpt.pos.x(5);   

th    = 0:4;

% SURGE EQUATION xd(th)

% velocity/acceleration at boundary
boundary = 'vel';   
%boundary  = 'acs'

start = 10;
final = 190;

y = [start x1  x2 x2 0 0  x3 x3 0 0  x4 x4 0 0  x5 final ]';

if boundary == 'vel',
    c1 = pva(th(1),'v');
    c2 = pva(th(5),'v');
elseif boundary == 'acs',
    c1 = pva(th(1),'a');
    c2 = pva(th(5),'a') 
end

O = zeros(1,4);

A = [ c1                O O O  
      pva(th(1),'p')    O O O                 
%
      pva(th(2),'p')               O  O O  
                   O  pva(th(2),'p')  O O   
      -pva(th(2),'v') pva(th(2),'v')  O O    
      -pva(th(2),'a') pva(th(2),'a')  O O 
%      
      O  pva(th(3),'p')               O  O 
      O              O   pva(th(3),'p')  O  
      O -pva(th(3),'v')  pva(th(3),'v')  O 
      O -pva(th(3),'a')  pva(th(3),'a')  O         
% 
      O  O pva(th(4),'p')                O   
      O  O              O   pva(th(4),'p')   
      O  O -pva(th(4),'v')  pva(th(4),'v')   
      O  O -pva(th(4),'a')  pva(th(4),'a')        
%
      O O O           pva(th(5),'p')   
      O O O           c2              ]; 

x = inv(A)*y;

Acoeff(1,:) = x(1:4)';  
Acoeff(2,:) = x(5:8)';  
Acoeff(3,:) = x(9:12)'; 
Acoeff(4,:) = x(13:16)'; 

% SWAY EQUATION yd(th)

y1     = wpt.pos.y(1);
y2     = wpt.pos.y(2);    
y3     = wpt.pos.y(3);   
y4     = wpt.pos.y(4);   
y5     = wpt.pos.y(5);   

th    = 0:4;

% velocity/acceleration at boundary
boundary = 'vel';
%boundary  = 'acs'

start = 0;
final = 0;

y = [start y1  y2 y2 0 0  y3 y3 0 0  y4 y4 0 0  y5 final ]';

if boundary == 'vel',
    c1 = pva(th(1),'v');
    c2 = pva(th(5),'v');
elseif boundary == 'acs',
    c1 = pva(th(1),'a');
    c2 = pva(th(5),'a') 
end

O = zeros(1,4);

A = [ c1                O O O  
      pva(th(1),'p')    O O O                 
%
      pva(th(2),'p')               O  O O  
                   O  pva(th(2),'p')  O O   
      -pva(th(2),'v') pva(th(2),'v')  O O    
      -pva(th(2),'a') pva(th(2),'a')  O O 
%      
      O  pva(th(3),'p')               O  O 
      O              O   pva(th(3),'p')  O  
      O -pva(th(3),'v')  pva(th(3),'v')  O 
      O -pva(th(3),'a')  pva(th(3),'a')  O         
% 
      O  O pva(th(4),'p')                O   
      O  O              O   pva(th(4),'p')   
      O  O -pva(th(4),'v')  pva(th(4),'v')   
      O  O -pva(th(4),'a')  pva(th(4),'a')        
%
      O O O           pva(th(5),'p')   
      O O O           c2              ]; 

x = inv(A)*y;

Bcoeff(1,:) = x(1:4)';  
Bcoeff(2,:) = x(5:8)';  
Bcoeff(3,:) = x(9:12)'; 
Bcoeff(4,:) = x(13:16)'; 

% GRAPHICS

figure(1)

subplot(321)
plot(0:4,wpt.pos.x(1:5),'ro')
hold on
for i= 1:4,
    th = linspace(i-1,i,101);  
    plot(th,polyval(Acoeff(i,:),th),'linewidth',2);
end
grid
title('x_d(\theta)')
hold off

subplot(325)
plot(0,0)
hold on
for i= 1:4,    
    th = linspace(i-1,i,101); 
    DDxth = polyval([0 0 6*Acoeff(i,1) 2*Acoeff(i,2) ],th);
    plot(th,DDxth,'linewidth',2);
end
hold off
title('x_d''''(\theta)')
grid 

subplot(323)

plot(0,0)
hold on
for i= 1:4,    
    th = linspace(i-1,i,101); 
    Dxth = polyval([0 3*Acoeff(i,1) 2*Acoeff(i,2) Acoeff(i,3)],th);
    plot(th,Dxth,'linewidth',2);
end
hold off
title('x_d''(\theta)')
grid 


subplot(322)
plot(0:4,wpt.pos.y(1:5),'ro')
hold on
for i= 1:4,
    th = linspace(i-1,i,101);  
    plot(th,polyval(Bcoeff(i,:),th),'linewidth',2);
end
grid
title('y_d(\theta)')
hold off

subplot(326)
plot(0,0)
hold on
for i= 1:4,    
    th = linspace(i-1,i,101); 
    DDyth = polyval([0 0 6*Bcoeff(i,1) 2*Bcoeff(i,2) ],th);
    plot(th,DDyth,'linewidth',2);
end
hold off
title('y_d''''(\theta)')
grid 

subplot(324)

plot(0,0)
hold on
for i= 1:4,    
    th = linspace(i-1,i,101); 
    Dyth = polyval([0 3*Bcoeff(i,1) 2*Bcoeff(i,2) Bcoeff(i,3)],th);
    plot(th,Dyth,'linewidth',2);
end
hold off
title('y_d''(\theta)')
grid 

figure(2)
subplot(111)
plot(wpt.pos.y(1:5),wpt.pos.x(1:5),'ro')
hold on
for i= 1:4,    
    th = linspace(i-1,i,101); 
    xth = polyval(Acoeff(i,:),th);
    yth = polyval(Bcoeff(i,:),th);
    plot(yth,xth,'linewidth',2);
end
hold off
xlabel('y_d(\theta)')
xlabel('x_d(\theta)')
axis('equal')
grid 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GUIDANCE CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear simdata

T    = 10; % speed time constant
Uref = 5;  % desired speed waypoint 1

U = 0;     % initial speed waypoint 0
th = 0;    % initial th-value waypoint 0

h = 0.1;   % sampling time
i = 0;     % counter

while th<1,
    
    i = i+1;
    
    % cubic polynominal between waypoints 0 and 1
    xth = polyval([0 3*Acoeff(1,1) 2*Acoeff(1,2) Acoeff(1,3)],th);
    yth = polyval([0 3*Bcoeff(1,1) 2*Bcoeff(1,2) Bcoeff(1,3)],th);
    
    th = th + h*(U / sqrt(xth^2+yth^2) );    % theta dynamics
    U  = U  + h*(-U + Uref)/T;               % speed dynamics
    
    % store solution of ODEs: simdata = [time, U, th]
    simdata(i,:) = [(i-1)*h,U,th];        
    
end

% graphics

t = simdata(:,1);
U = simdata(:,2);
th = simdata(:,3);

figure(3)

subplot(211)
plot(t,ones(size(t))*Uref)
hold on
plot(t,U,'linewidth',2)
hold off
grid
title('Speed U_d(t) as a function of time t')

subplot(212)
plot(t,ones(size(t)))
hold on
plot(t,th,'linewidth',2),
hold off
grid
title('Path variable \theta(t) as a function of time t')