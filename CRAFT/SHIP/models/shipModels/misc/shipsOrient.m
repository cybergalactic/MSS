% Vessel CG track and shape plotter
% 
% B     = vessel breadth (m)
% L     = vessel length  (m)
% 
% B & L should be loaded from file with vessel's particulars
% Other data is captured from Simulink model
%
% Vessel shape is determined as five consequently connected points 
% determinned from B & L.
%
% K     = vessel course             (rad)
% XY    = vessel CG track [x,y]
% XY1   = vessel CG track, taken for shape plotting period 
% (specified in Simulink model)
%
% Author:    Aleksandr D. Pipchenko
% Date:      20th June 2007
% Revisions: 6 March 2011 (Aleksandr D. Pipchenko): added description 

XX=[0 B/2 B/2 -B/2 -B/2 0];YY=[L/2 L/3.5 -L/2 -L/2 L/3.5 L/2];
Ao=atan(XX./YY);
R = YY./cos(Ao);
XC =[-B/2 -B/2];YC = [-60 40];Ac=atan(XC./YC);Rc = YC./cos(Ac);

for i = 1:numel(XY1(:,1))
    for j = 1:numel(XX)
        A(j,i)=Ao(j)+K(i);
        XN(j,i)=XY1(i,1)+R(j)*sin(A(j,i));YN(j,i)=XY1(i,2)+R(j)*cos(A(j,i));
    end
end

close all; figure(10);
plot (XN,YN,'b',XY(:,1),XY(:,2),'g');grid on;axis equal;
clear XN; clear YN; clear XY1;