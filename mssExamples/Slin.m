function Pyy = Slin(lambda,w)
% Pyy = Slin(lambda,w) 2nd-order linear power spectracl density (PSD) function
%
%   w       = wave spectrum frequency (rad/s)
%   lambda  = relative damping factor
%
%   See ExLinspec.m
%
% Author:   Thor I. Fossen
% Date:     15th August 2001
% Revisions: 

global sigma wo  

Pyy = 4*(lambda*wo*sigma)^2*w.^2 ./ ( (wo^2-w.^2).^2 + 4*(lambda*wo.*w).^2 );

