echo off
% WAVEDEMO  Wave spectrum demonstration
%
% Author:     Thor I. Fossen
% Date:       2001-08-14
% Revisions:  

echo on; close all

Hs   = 10;       % significant wave height
wmax = 3.0;      % maximum wave frequency
To   = 5.0;      % peak period of the wvae spectrum
wo   = 2*pi/To;  % peak frequency 

%  SpecType and Par =[p1,p2,p3,...pn]:
%   SpecType =1  Bretschneither (p1=A,p2=B)
%   SpecType =2  Pierson-Moskowitz   (p1=Vwind20) 
%   SpecType =3, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T0)
%   SpecType =4, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T1)
%   SpecType =5, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=Tz)
%   SpecType =6, JONSWAP (p1=Vwind10,p2=Fetch)
%   SpecType =7, JONSWAP (p1=Hs,p2=w0,p3=gamma)
%   SpecType =8, Torsethaugen (p1=Hs,p2=w0) 
  
w = (0:0.025:2)';
S1 = wavespec(3,[Hs,To],w,0);       % modified Pierson-Moskowitz spectrum
S2 = wavespec(7,[Hs,wo,3.3],w,0);   % JONSWAP gamma = 3.3
S3 = wavespec(7,[Hs,wo,2.0],w,0);   % JONSWAP gamma = 2.0
S4 = wavespec(8,[Hs,wo],w,0);       % Torsethaugen

plot(w,S1/(To*Hs^2),'k*-',w,S2/(To*Hs^2),'k-',w,S3/(To*Hs^2),'k-.',w,S4/(To*Hs^2),'ko-')
legend('Modified Pierson-Moskowitz','JONSWAP for \gamma = 3.3',...
    'JONSWAP for \gamma = 2.0','Torsethaugen','Location','Northeast')
xlabel('\omega (rad/s)')
title('S(\omega)/(H_s^2 T_0)')
grid
axis([0 2 0 0.04])

set(gca,'FontSize',12)

pause % Strike any key to return to MAIN MENU

close(1); echo off
disp('End of demo')

