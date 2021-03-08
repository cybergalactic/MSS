disp('...running cylinderModel.m')
clear all
load cylinder
load cylinderABC

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
    T(i) = 2*pi/w_n;
end

disp('Natural periods:')
disp(T)



function F = natfreq(x,m,k,w)
   F = x - sqrt(k/interp1(w,m,x));
end



