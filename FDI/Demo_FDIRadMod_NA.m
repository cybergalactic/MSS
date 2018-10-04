
%Demo: example use of frequency-domain identification of radiation models
%of marine strucutres based on hydrodynamic data without using the
%infinite-frequency added mass.
%
% Created by Tristan Perez (tristan.perez@ntnu.no)
% Date 2009/9/1, Trondheim, Norway.
% Revision:

clear all

%Load data of the fpso (structure vessel, 6DOF)
load fpso.mat

%Extract the data from the vessel structure

Dof = [3,3]; %Use coupling 3-3 heave-heave
Nf = length(vessel.freqs);
W=vessel.freqs(1:Nf-1)';
Ainf=vessel.A(Dof(1),Dof(2),Nf);

A = reshape(vessel.A(Dof(1),Dof(2),1:Nf-1),1,length(W))';
B = reshape(vessel.B(Dof(1),Dof(2),1:Nf-1),1,length(W))';

%Define the structure with identification algorithm options
FDIopt.OrdMax     = 20;
FDIopt.AinfFlag   = 0;
FDIopt.Method     = 2;
FDIopt.Iterations = 20;
FDIopt.LogLin     = 1;
FDIopt.wsFactor   = 0.1;  
FDIopt.wminFactor = 0.1;
FDIopt.wmaxFactor = 5;

%call idenfication routine
[Krad,Ainf_hat]=FDIRadMod(W,A,0,B,FDIopt,Dof);


