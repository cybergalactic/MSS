function [We]=ww2we(chi,U,Ww)
% Function to transform from a vector wave frequancy to 
% encounter frequencies
%
%Use: [We]=ww2we(chi,U,W)
% 
% We   vecrtor of encounter fequency values [rad/sec] 
% chi  encounter angle [rad] 0-following seas, pi-head seas
% U    forward speed [m/sec]
% W    vecrtor of wave fequancy values [rad/sec] 
%
%Reference: A.R.M.J LLoyd "Seakeeping: Ship Behaviour in Rough Wheather."
%John Wiley & Sons, 1989.
%
% Created by: Tristan Perez in 2001  
% Last mod. by: Tristan Perez 
% Date: 9 March 2005

w=[];
for k=1:length(Ww)
      we_aux=(Ww(k)-Ww(k)^2*U*cos(chi)/9.81);
   w=[w;we_aux];
end
We=w;