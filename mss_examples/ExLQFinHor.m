% ExLQFinHor LQ finite time-horizon tracking of mass-damper-spring system
% Author:    Roger Skjetne
% Date:      28th Januray 2001
% Revisions: 15th June 2001, Thor I. Fossen - minor changes of notation 

format long;

% Time horizon
T = 10;
Ts = 0.01;
N = round(T/Ts);

% Plant:
A = [0 1; -1 -2]; B = [0; 1]; E = [0; 1]; C = [1 0];
plantC = ss(A,[B E],C,zeros(1,length([B E])));  % plantC is of type structure
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
for k=1:N+1,
    Xref(:,k) = xr;                    % Storing the values in matrix Xref
    xr = RefGenD.a*xr + RefGenD.b*sin((k-1)*Ts);    % Updating
end

% Backwards simulation. Let Xc = [xd' Pvec' h1' h2']'
PT = C.'*Qf*C; P = PT;
h1T = -C.'*Qf*C*xr0;  h1 = h1T;
h2T = zeros(size(B)); h2 = h2T;

H1 = zeros(length(h1T),N+1);
H2 = zeros(length(h2T),N+1);
Ps = struct('Mat',[]);

for k=1:N+1,
    H1(:,(N+1)-(k-1)) = h1;             % Storing the values in matrix H1
    H2(:,(N+1)-(k-1)) = h2;             % Storing the values in matrix H2
    Ps(k) = struct('Mat',P);            % Storing the values in structure Ps
    
    K = A - B*inv(R)*B.'*P;
    h1dot = -K.'*h1 + Qt*Xref(:,(N+1)-(k-1));
    h2dot = -K.'*h2 - P*E*cos((N-(k-1))*Ts);
    Pdot = -P*A - A.'*P + P*B*inv(R)*B.'*P - Qt;
    
    h1 = h1 - Ts*h1dot;                 % New previous value for h1
    h2 = h2 - Ts*h2dot;                 % New previous value for h2
    P = P - Ts*Pdot;                    % New previous value for P
end


for k=1:N+1,
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
    
    u(k) = -inv(R)*B.'*(Ps((N+1)-(k-1)).Mat*x + H1(:,k) + H2(:,k));
    if abs(u(k)) > Umax
        u(k) = sign(u(k))*Umax;
    end
    x = plantD.a*x + plantD.b*[u(k); cos((k-1)*Ts)];    % Updating
end

% Plotting the results
figure(1)
k = 0:N;
subplot(211);
plot(k*Ts, X, k*Ts, Xref, '-.');
legend('x_1(t)','x_2(t)','x_{r1}(t)','x_{r2}(t)');
title('Output response to input r(t) = sin(t)');
xlabel('Time (sec.)'); 

subplot(212);
plot(k*Ts, u);
legend('u(t)'); axis([0 T -Umax-1 Umax+1]);
xlabel('Time (sec.)'); ylabel('Control input');

% Plotting the results:
figure(2)
k = 0:N;
subplot(311);
plot(k*Ts, p11,':', k*Ts, p12,'-.', k*Ts, p22,'--');
legend('p_{11}(t)','p_{12}(t)','p_{22}(t)');
title('Dynamic behavior of P(t)');
xlabel('Time (sec.)');

subplot(312);
plot(k*Ts, H1(1,:),':', k*Ts, H1(2,:),'-.');
legend('h_{11}(t)', 'h_{12}(t)'); 
title('Dynamic behavior of h_1(t)');
xlabel('Time (sec.)')

subplot(313);
plot(k*Ts, H2(1,:),':', k*Ts, H2(2,:),'-.');
legend('h_{21}(t)', 'h_{22}(t)'); 
title('Dynamic behavior of h_2(t)');
xlabel('Time (sec.)')
