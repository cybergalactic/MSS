disp('running model.m...')
clear all
load cylinder

% mooring
C_mooring = diag([50e3 50e3 0 0 0 397e3]);

% constants
mass = vessel.main.m;
g    = vessel.main.g;

% natural periods in CO
for i  = 1:6
    k = vessel.C(i,i,end)+C_mooring(i,i);
    m = vessel.MRB(i,i)+reshape(vessel.A(i,i,:),1,length(vessel.freqs));
    w = vessel.freqs;
    [w_n,F_n,flag] = fsolve(@(x) natfreq(x,m,k,w),0.5,optimset('Display','off','TolFun',1e-10));  
    w_o(i) = w_n;
    T(i) = 2*pi/w_n;
end

disp('Natural periods:')
disp(T)

% Simulink matrices
idx  = 25;
w_design = vessel.freqs(25)            % frequency used for data

MRB  = vessel.MRB;
MA   = vessel.A(:,:,idx);
M    = MRB + MA;
Minv = inv(M);

D    = vessel.B(:,:,idx);
G    = vessel.C(:,:,end);


% add linear viscous damping in percentage (relative damping factors)
D(1,1) = D(1,1) + 2*0.05*w_o(1)*M(1,1);
D(2,2) = D(2,2) + 2*0.05*w_o(2)*M(2,2);
D(6,6) = D(6,6) + 2*0.10*w_o(6)*M(6,6);

D(3,3) = D(3,3) + 2*0.05*w_o(3)*M(3,3);
D(4,4) = D(4,4) + 2*0.05*w_o(4)*M(4,4);
D(5,5) = D(5,5) + 2*0.05*w_o(5)*M(5,5);

