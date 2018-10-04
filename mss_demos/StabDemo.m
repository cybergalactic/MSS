% STABDEMO     Straight-line, directional and positional motion stability
% Author:      Thor I. Fossen
% Date:        2001-07-08
% Revisions:   

N = 6000;  % number of samples
h = 0.1;   % sampling time

% cargo ship
K = 0.185;  T = 107.3;  uo = 5.0;

No = menu('Choose maneuver','Straight-line stability',...
          'Directional stability (critical damped)',...
          'Directional stability (underdamped)',...
          'Positional motion stability',...
          'Exit');
    
% main loop
r=0; psi=0; x=0; y=0; delta = 0; z=0;
xout = zeros(N+1,5);

for i=1:N+1,
    xout(i,:) = [r psi x y delta];
     
    if No==1,
        delta = 0;
    elseif (No==2 | No==3),
        if No==2, zeta = 1;   wn = 3/T;  end
        if No==3, zeta = 0.1; wn = 3/T;  end        
        Kp = (T/K)*wn*wn;              % PD-control
        Kd = (T/K)*(2*zeta*wn+1/T);
        delta = -Kp*(psi)-Kd*r;
    elseif No==4,
        zeta = 1; wn = 7/T;  
        Kp = (T/K)*wn*wn; 
        Kd = (T/K)*(2*zeta*wn+1/T);
        Ki = Kp*wn/10;
        delta = -Kp*(psi)-Kd*r-Ki*z;   % PID-control
        z = z + h*psi;
    else
        break
    end

    w = 0;
    if i==1000, w = 1; end   % impulse w(t) at time t = 100 (s)
    
    r   = r + h*(-r + K*delta + w)/T;  % Euler integration
    psi = psi + h*r;
    x = x + h*uo*cos(psi);
    y = y + h*uo*sin(psi);
end

% plots
if No~=5,  
    t = h*(0:N)';
    r=(180/pi)*xout(:,1); psi=(180/pi)*xout(:,2); x=xout(:,3); y=xout(:,4); delta=(180/pi)*xout(:,5);

    clf;figure(1),figure(gcf)
    subplot(211);plot(x,y,'b','linewidth',2),hold on 
    plot(x,0*zeros(length(x),1),'r','linewidth',2),hold off
    grid
    if No==1, title('XY-Plot: Straight-line stability'); 
    elseif No==2, title('XY-Plot: Directional stability (critical damped)');
    elseif No==3, title('XY-Plot: Directional stability (underdamped)');
    elseif No==4, title('XY-Plot: Positional motion stability');
    end
    subplot(223);plot(t,r,'linewidth',2),grid,title('r (deg/s)'),xlabel('sec')
    subplot(224);plot(t,psi,'linewidth',2),grid,title('\psi (deg)'),xlabel('sec')

    StabDemo
else
    close all
end