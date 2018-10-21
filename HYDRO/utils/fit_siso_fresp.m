function [Num,Den,Fresp_hat,stable]=fit_siso_fresp(Omega,Fresp,Weight,ord_num,ord_den,method,iter)
% Function to fit a continuous SISO TF to 
% a frequency response data based on invfreqs.m 
%
% The non parametric
% data H(jw) is used to estimate the parameters of a parametric model
% H(s) = P(s)/Q(s), and the parameters of P(s) and Q(s) are estimated
% from a LS problem:
%
% Theta_opt = argmin Sum_i |Kw(jwi)-P(Theta,jwi)/Q(Theta,jwi)|^2
%
%
% Use:
% [Num,Den,Fresp_hat]=fit_siso_fresp(Omega,Fresp,Weight,ord_num,ord_den,method,iter)
%
% inputs:
% Omega -  vector of freqs rad/sec
% Fresp -  complex freq response of the system to approximate
% Weight-  vector of wightings for LS
% ord_num- order of numerator
% ord_den- order of denominator
% iter   - number of iterations
% method - =1 Levy (quasi linearization) Theta_opt = argmin Sum_i |Q(Theta,jwi) Kw(jwi)-P(Theta,jwi)|^2   
%          =2 iterative Levy with Den from previous iterations as previous weight
%          =3 NL LS based on Gauss-Newton
%
% The method 2 is faster than 3 and gives the same results.
%s
% Author: Tristan Perez 
% Date:   2007.8.4
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

if method==1,
    
    [Num,Den] = invfreqs(Fresp,Omega,ord_num,ord_den,Weight);
    Fresp_hat = freqs(Num,Den,Omega);
    stable    = MyIsStable(Den);

elseif method==2,
    
    Weight = ones(size(Omega));
    for k=1:iter,
        [Num,Den] = invfreqs(Fresp,Omega,ord_num,ord_den,Weight);
        for m=1:length(Omega);
            Weight(m)=1/abs(polyval(Den,complex(0,Omega(m))))^2;
        end
    end
    Fresp_hat=freqs(Num,Den,Omega);
    stable = MyIsStable(Den);
    
elseif method==3,
    
    [Num,Den] = invfreqs(Fresp,Omega,ord_num,ord_den,Weight,iter);
    Fresp_hat = freqs(Num,Den,Omega);
    stable    = MyIsStable(Den);
    
else
    disp('Error: Method is 1,2, or 3')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = MyIsStable(p)
% This function checks if a polynomial has roots in the 
% closed right half plane, and if so, the unstable roots are 
% reflected about the imaginary axis.
ut=1;
r=roots(p);
for k=1:length(r)
    if  real(r(k)) >= 0,
        ut=0;
    end
end
out=ut;