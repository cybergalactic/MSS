% ExLQFinHor LQ finite time-horizon tracking of mass-damper-spring system
% Author:    Roger Skjetne
% Date:      28 Januray 2001
% Revisions: 15 June 2001, T. I. Fossen - minor changes of notation 
%             9 June 2020, T. I. Fossen - removed disturbance FF
format long;

% Time horizon
T = 10;
Ts = 0.01;
N = round(T/Ts);

% Plant:
A = [0 1; -1 -2]; B = [0; 1]; E = [0; 0]; C = [1 0];
plantC = ss(A,[B E],C,zeros(1,length([B E])));  % plant C is of type structure
[plantD, Mp] = c2d(plantC, Ts, 'foh');          % discretization
x0 = [1; 2];                                    % initial states

Umax = 20;                                      % Saturation level: -Umax < u < +Umax.

% Reference system:
Ad = [0 1; -1 -1]; Bd = [0; 1];
RefGenC = ss(Ad,Bd,C,0);
[RefGenD, Mr] = c2d(RefGenC, Ts, 'foh');    % RefGenD is of type structure
xr0 = [0; 0];                               % initial states

% Weight Matrices:
Q  = 1000; Qt = C.'*Q*C; R  = 5; Qf = 0;

% Simulate the reference system forward to obtain all 
% future knowledge of this system up to time T.     r(t) = sin(t):
xr = xr0;
Xref = zeros(length(xr0),N+1);
for k=1:N+1
    Xref(:,k) = xr;                    % Storing the values in matrix Xref
    xr = RefGenD.a*xr + RefGenD.b*sin((k-1)*Ts);    % Updating
end

% Backwards simulation. Let Xc = [xd' Pvec' h']'
PT = C.'*Qf*C; P = PT;
hT = -C.'*Qf*C*xr0;  h = hT;

H = zeros(length(hT),N+1);
Ps = struct('Mat',[]);

for k=1:N+1
    H(:,(N+1)-(k-1)) = h;               % Storing the values in matrix H
    Ps(k) = struct('Mat',P);            % Storing the values in structure Ps
    
    K = A - B*inv(R)*B.'*P;
    hdot = -K.'*h + Qt*Xref(:,(N+1)-(k-1));
    Pdot = -P*A - A.'*P + P*B*inv(R)*B.'*P - Qt;
    
    h = h - Ts*hdot;                    % New previous value for h
    P = P - Ts*Pdot;                    % New previous value for P
end


for k=1:N+1
    p11(k) = Ps((N+1)-(k-1)).Mat(1);
    p12(k) = Ps((N+1)-(k-1)).Mat(2);
    p22(k) = Ps((N+1)-(k-1)).Mat(4);
end

% Final Simulation of the total system: Total state X = [x' Xc']'
% Let plant Initial Condition be x0 = [1 2]'
x = x0; u = 0;
X = zeros(length(x0),N+1);

for k = 1:N+1
    X(:,k) = x;                         % Storing the values in matrix X
    
    u(k) = -inv(R)*B.'*(Ps((N+1)-(k-1)).Mat*x + H(:,k));
    if abs(u(k)) > Umax
        u(k) = sign(u(k))*Umax;
    end
    x = plantD.a*x + plantD.b*[u(k); cos((k-1)*Ts)];    % Updating
end

% Plotting the results
figure(1)
k = 0:N;
subplot(411);
plot(k*Ts, X, k*Ts, Xref, '-.','linewidth',2);
legend('x_1','x_2','x_{r1}','x_{r2}');
title('Output response to input r(t) = sin(t)');
xlabel('time (s)'); 

subplot(412);
plot(k*Ts, u,'linewidth',2);
legend('u'); axis([0 T -Umax-1 Umax+1]);
xlabel('time (s)'); ylabel('Control input');

% Plotting the results:
subplot(413);
plot(k*Ts, p11,'k-', k*Ts, p12,'-.', k*Ts, p22,'k--','linewidth',2);
legend('p_{11}','p_{12}','p_{22}');
title('Dynamic behavior of P(t)');
xlabel('time (s)');

subplot(414);
plot(k*Ts, H(1,:),'-.', k*Ts, H(2,:),'-','linewidth',2);
legend('h_{11}', 'h_{12}'); 
title('Dynamic behavior of h(t)');
xlabel('time (s)')

