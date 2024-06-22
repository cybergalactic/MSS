echo off
% WAVEDEMO: MSS wave spectrum demonstration.
%
% Author:     Thor I. Fossen
% Date:       2001-08-14
% Revisions:  

echo on; close all

Hs   = 10;       % Significant wave height
wmax = 3.0;      % Maximum wave frequency
To   = 5.0;      % Peak period of the wvae spectrum
wo   = 2*pi/To;  % Peak frequency 

%  SpecType and Par = [p1, p2, p3,..., pn]:
%   SpecType = 1 : Bretschneither (p1 = A, p2 = B)
%   SpecType = 2 : Pierson-Moskowitz (p1 = Vwind20) 
%   SpecType = 3 : ITTC-Modified Pierson-Moskowitz (p1 = Hs, p2 = T0)
%   SpecType = 4 : ITTC-Modified Pierson-Moskowitz (p1 = Hs, p2 = T1)
%   SpecType = 5 : ITTC-Modified Pierson-Moskowitz (p1 = Hs, p2 = Tz)
%   SpecType = 6 : JONSWAP (p1 = Vwind10, p2 = Fetch)
%   SpecType = 7 : JONSWAP (p1 = Hs, p2 = w0, p3 = gamma)
%   SpecType = 8 : Torsethaugen (p1 = Hs, p2 = w0) 
  
w = (0:0.025:wmax)';
S1 = wavespec(3,[Hs,To],w,0);       % Modified Pierson-Moskowitz spectrum
S2 = wavespec(7,[Hs,wo,3.3],w,0);   % JONSWAP with gamma = 3.3
S3 = wavespec(7,[Hs,wo,2.0],w,0);   % JONSWAP with gamma = 2.0
S4 = wavespec(8,[Hs,wo],w,0);       % Torsethaugen

plot(w,S1/(To*Hs^2),'*-',w,S2/(To*Hs^2),'-',w,S3/(To*Hs^2),'-.', ...
    w,S4/(To*Hs^2),'o-')
legend('Modified Pierson-Moskowitz','JONSWAP for \gamma = 3.3',...
    'JONSWAP for \gamma = 2.0','Torsethaugen','Location','Northeast')
xlabel('\omega (rad/s)')
title('S(\omega)/(H_s^2 T_0)'); grid; echo off

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)



