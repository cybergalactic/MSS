% plotAterms plots the CA terms (speed dependent terms)
% _________________________________________________________________________
%
% MSS HYDRO is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

figure(6)

% SS11 = ss(AAr2(:,:,3,3),BBr2(:,:,3,3),CCr2(:,:,3,3),DDr(:,:,3,3));
% SS22 = ss(AAr2(:,:,3,5),BBr2(:,:,3,5),CCr2(:,:,3,5),DDr(:,:,3,5));
% SS33 = ss(AAr2(:,:,2,2),BBr2(:,:,2,2),CCr2(:,:,2,2),DDr(:,:,2,2));
% SS44 = ss(AAr2(:,:,2,4),BBr2(:,:,2,4),CCr2(:,:,2,4),DDr(:,:,2,4));
% SS55 = ss(AAr2(:,:,1,1),BBr2(:,:,1,1),CCr2(:,:,1,1),DDr(:,:,1,1));
% SS66 = ss(AAr2(:,:,2,6),BBr2(:,:,2,6),CCr2(:,:,2,6),DDr(:,:,2,6));

ARplot = reshape(A_retard(3,3,:),1,length(t));
subplot(3,2,1)
plot(t,ARplot,'b-',t,imp(:,1),'r'),grid
ylabel('A_{33}'),xlabel('time (s)')

ARplot = reshape(A_retard(3,5,:),1,length(t));
subplot(3,2,2)
plot(t,ARplot,'b-',t,imp(:,2),'r'),grid
ylabel('A_{35}'),xlabel('time (s)')

ARplot = reshape(A_retard(2,2,:),1,length(t));
subplot(3,2,3)
plot(t,ARplot,'b-',t,imp(:,3),'r'),grid
ylabel('A_{22}'),xlabel('time (s)')

ARplot = reshape(A_retard(2,4,:),1,length(t));
subplot(3,2,4)
plot(t,ARplot,'b-',t,imp(:,4),'r'),grid
ylabel('A_{24}'),xlabel('time (s)')

ARplot = reshape(A_retard(1,1,:),1,length(t));
subplot(3,2,5)
plot(t,ARplot,'b-',t,imp(:,5),'r'),grid
ylabel('A_{11}'),xlabel('time (s)')

ARplot = reshape(A_retard(2,6,:),1,length(t));
subplot(3,2,6)
plot(t,ARplot,'b-',t,imp(:,6),'r'),grid
ylabel('A_{26}'),xlabel('time (s)')