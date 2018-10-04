function msi = HMmsi(a_z,w_e)
% msi = HMMSI(a_z,w_e) computes the Motion Sickness Incidence using the method 
%       of O'Hanlon and McCauley (1974).
%
% a_z: mean of absolute vertical acceleration (m/s^2), i.e. a_z = mean(abs(a_measured(:,1)))
% w_e: vector of enconter frequency (rad/s), see encounter.m.
% msi: the percentage of persons that become seasick during a 2 hours sail
%
% Refs.  - A. R. J. M. Lloyd (1989). Seakeeping Behaviour in Rough Water. Ellis Horwoowd Ltd.
%        - E. V. Lewis (Ed.) (1989). Principles of Naval Architecture. Vol III Motions in
%             Waves and Controllability, 2nd. ed., SNAME. 
%        - J.F. O'Hanlon and M. E. McCauley (1974). Motion Sickness Incidence as a Function of
%             Vertical Sinusoidal Motion. Aerospace Medicine AM-45(4):366-369.
%
% Author:    J. M. de la Cruz
% Date:      28th March 2000
% Revisions: 5th November 2001, Thor I. Fossen - minor changes of notation/documentation 

mu_MSI = -0.819+2.32*(log10(w_e)).^2;
I      = ( -log10(abs(a_z)/9.8 )+mu_MSI )/0.4;

% Divide by sqrt(2) to compensate the difference in the definition of the 
% Matlab function erf.m and erf() used in Lloyd (1989), page 335.
I = I/sqrt(2); 

for i=1:length(w_e),
    if I(i) >=0,
        msi(i)=0.5-0.5*abs(erf(I(i)));
    else
        msi(i)=0.5+0.5*abs(erf(I(i)));
    end
end

msi=100*msi;