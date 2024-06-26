% exRefMod   2nd-order reference model with nonlinear damping and velocity 
% saturation.
%
% Author:    Thor I. Fossen
% Date:      3rd November 2001
% Revisions: 

zeta = 1;  % Relative damping ratio
w = 1;     % Natural frequency
delta = 1; % Nonlinear damping coefficients
v_max = 1; % Max velocity
h = 0.1;   % Sampling time
N = 200;   % Number of samples

t  = 0:h:h*N;

% Linear mass-damper-spring system
r1 = 10 * ones(max(size(t)),1);
[A,B,C,D] = ord2(w,zeta); 
[x1,y1]   = lsim(A,B,C,D,r1,t); 

r2 = 10*r1;
[A,B,C,D] = ord2(w,zeta); 
[x2,y2]   = lsim(A,B,C,D,r2,t); 

% Nonlinear damping
x = 0;
v = 0;
r2 = 10;
y2 = zeros(N+1,2);
for i=1:N+1
   y2(i,:) = [x v];   
   x_dot = v;
   v_dot = w^2 * (r2 - x) - 2 * zeta * w * v - delta *abs(v) * v;
   v = v + h * v_dot;
   x = x + h * x_dot;
end

% Velocity saturation
x = 0;
v = 0;
r3 = 10;
y3 = zeros(N+1,2);

for i=1:N+1

   y3(i,:) = [x v];   

   v_dot = w^2 * (r3 -x) - 2 * zeta * w * v;
   x_dot = v;

   v = v + h * v_dot;
   if abs(v) > v_max           % Velocity saturation
      v = sign(v) * v_max;
   end
   x = x + h * x_dot;

end

%% Plots
figure(gcf)
subplot(211); plot(t,y1(:,1))
hold on; plot(t,y2(:,1),'--k',t,y3(:,1),'r'); hold off; grid
title('2nd-order mass-damper-spring reference model')
legend('Linear damping','Nonlinear damping','Velocity saturation')
subplot(212); plot(t,y1(:,2));
hold on; plot(t,y2(:,2),'--k',t,y3(:,2),'r'); hold off; grid

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)


