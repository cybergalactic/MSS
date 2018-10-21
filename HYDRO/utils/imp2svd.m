function [SVD]=imp2svd(Y)
% Y - is verctor containg the samples of the impulse
%       response of a SISO system.
%
% Author: Tristan Perez 
% Date:   2007.8.4
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).

[mm,nn] = size(Y);
if nn==1,
  Y=Y';
  [mm,nn] = size(Y);
end
nt=nn-1;
h=hankel(Y(2:nn));
[u,SVDM,v] = svd(h);
SVD=SVDM*ones(max(size(SVDM)),1);