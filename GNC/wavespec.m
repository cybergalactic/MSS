function S = wavespec(SpecType,Par,W,PlotFlag)
%Function to evaluate different type of wave spectra.
%
%Use:  S = wavespec(SpecType,Par,W,PlotFlag)
%
%Input:
% SpecType	- Spectrum type
% Par		- Spectrum parameters
% W			- Column vector of wave frequencies [rad/s]
% PlotFlag	- 1 to plot the spectrum, 0 for no plot
%
%Output:
% S			- Column vector of wave spectrum values [m^2 s]; evaluated at W(k) 
%     
%SpecType and Par =[p1,p2,p3,...pn]:
% SpecType =1  Bretschneitder (p1=A,p2=B)
% SpecType =2  Pierson-Moskowitz   (p1=Vwind20) 
% SpecType =3, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T0)
% SpecType =4, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T1)
% SpecType =5, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=Tz)
% SpecType =6, JONSWAP (p1=Vwind10,p2=Fetch)
% SpecType =7, JONSWAP (p1=Hs,p2=w0,p3=gamma)
% SpecType =8, Torsethaugen (p1=Hs,p2=w0) 
%
%Bretschneither:
% p1 = A
% p2 = B
% S(w)=A*w^(-5)*exp(-B/(w^4));
% Reference [6]
%
%Pierson-Moskowitz:
% p1 = Vwind20 --Average wind speed @20m above sea level [m/s]
% A=8.1e-3*9.81^2;  
% B=0.74*(9.81/Vwind20)^4;
% S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
% Reference [5,6]
%
%
%ITTC-Modified Pierson-Moskowitz (Hs,T0):
% p1 = Hs  --Significant wave height (Hs = 4 m0^1/2) [m]
% p2 = T0  --Modal period      (T0 = 2 pi /w0)       [s]
% A=487*Hs^2/T0^4;
% B=1949/T0^4;
% S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
% Reference [1]
%
%
%ITTC-Modified Pierson-Moskowitz (Hs,T1):
% p1 = Hs  --Significant wave height (Hs = 4 m0^1/2) [m]
% p2 = T1   --Average wave period (T1 = 2 pi m0/m1)     [s]
% A=173*Hs^2/T1^4;
% B=691/T1^4;
% S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
% Reference [1]
%
%
%ITTC-Modified Pierson-Moskowitz (Hs,Tz):
% p1 = Hs --Significant wave height (Hs = 4 m0^1/2)             [m]
% p2 = Tz --Average zero-crossing period (T = 2 pi (m0/m1)^1/2) [s]
% A=123*Hs^2/Tz^4;
% B=495/Tz^4;
% S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
% Reference [1]
%
%
%JONSWAP (Vwind10,Fetch):
% p1 = Vwind10 --wind speed @ 10m over sea surface [m/sec]
% p2 = fetch --distance to georaphical boundary [m]
% g=9.81;
% xtilde= g*fetch/(Vwind10^2);
% f0=3.5*(g/Vwind10)*xtilde^-0.33;
% w0=2*pi*f0;
% alpha=0.076*xtilde^-0.22;
% gamma =3.3;
% sigma=0.07  if w<w0, sigma=0.09 otherwise;
% S(w)=S1*S2 -  [m^2 sec]  
%with,
% S1=alpha*g^2*(W^-5)*exp(-(5/4)*(w0/w)^4);
% S2=gamma^(exp(-(w-w0)^2/(2*(sigma*w0)^2)));   
% Reference [2]
%
%
%JONSWAP (Hs,w0, gamma):
% p1 = Hs - Significant wave height (Hs = 4 sqrt(m0)) [m]
% p2 = w0 - Modal Freq. [rad/sec] (Recomended 1.25<w0*sqrt(Hc)<1.75)
% p3 = gamma - Peakedness factor (Recommended between 1 and 5; usually 3.3, set to zero to use DNV formula)
% alpha=0.2*Hc^2*w0^4/g^2;
% g=9.81 [m/s^2]
% sigma=0.07  if w<w0, sigma=0.09 otherwise;
% S(w)=S1*S2  [m^2 s]  
%with,
% S1=alpha*g^2*(W^-5)*exp(-(5/4)*(w0/w)^4);
% S2=gamma^(exp(-(w-w0)^2/(2*(sigma*w0)^2))); 
% Reference [3]
%
%
%Torsethaugen (Hs,w0):
% The Torsethaugen spectrum is an empirical two peaked spectrum for swell 
% and developing sea based on experimental data from the North Sea. For 
% small peak frequencies, i.e.  0 < wmax <= 0.6 only one peak in the 
% spectrum appears. Returns the spectral density function S of the 
% Torsethaugen spectrum for the frequencies in the vector  W [rad/s].
%
% p1 = Hs  --significant wave height [m] 
% p2 = w0  -- Modal (peak) frequency [rad/s]
% S  - vector of power spectral densities [m^2s]
% Reference [4]
%
%
% References: 
%
% [1] A.R.M.J LLoyd (1998) "Seakeeping: Ship Behaviour in Rough Wheather."
%     Published by A.R.M.J LLoyd, Gosport UK.ISBN 0-9532634-01
%
% [2] Ochi, M.K. (1998) "Ocean Waves, The stochastic Approach"
%     Cambridge Ocean Technology Series, Vol 6,
%     Cambridge University Press.
%
% [3] S�rensen, A.J. (2005) "Marine Cybernetics: Modelling and Control"
%     Lecture Notes for TMR4241 Marine Control Systems, NTNU, 2005.
%
% [4] K.Torsethaugen (1996): "Model for a Doubly Peaked Wave Spectra"
%      Sintef report no.: STF22 A96204 prepared for Norsk Hydro.
%
% [5] T.I. Fossen (2002) "Marine Control Systems" Marine Cybernetics. 
%
% [6] Lewis E.V. "Principles of Naval Architecture volume  III
%     Motions in Waves and Controllability." SNAME, 1989.
%
% Created by: Tristan Perez in 2005  
% Revisions: 2007-03-09 minor fixes �yvind Smogeli
%             Adapted for MSS V3.0 March 2005. Grouped all spectra in one
%             function. Added Spec_type=1,4,5,6.
%            2009-09-11 fixed scaling of JONSWAP sepctrum

S=[];
[m,n] = size(W);
if n>m
   disp('Error: W must be a column vector')
   return
end

switch SpecType,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
case 1, %Bretschneither
	A = Par(1);
    B = Par(2);
    for k=1:1:length(W),
        Sa=A*W(k)^(-5)*exp(-B/(W(k)^4));
        S=[S;Sa];
    end
    TitleStr='Bretschneither Spectrum';
    L1Str = ['A=',num2str(A),' [m^2 s^{-4}], ','B=',num2str(B),' [s^{-4}]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2,%Pierson-Moskowitz 
	Vwind20 = Par(1);
    A=8.1e-3*9.81^2;  
    B=0.74*(9.81/Vwind20)^4;
    for k=1:1:length(W),
        Sa=A*W(k)^(-5)*exp(-B/(W(k)^4));
        S=[S;Sa];
    end
    TitleStr='Pierson-Moskowitz Spectrum';
    L1Str = ['Vwind @20m ASL =',num2str(Vwind20),' [m/s]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 3,%ITTC-Modified Pierson-Moskowitz (Hs,T0)
	Hs = Par(1);
    T0 = Par(2); 
    A=487*Hs^2/T0^4;
    B=1949/T0^4;  
    for k=1:1:length(W),
        Sa=A*W(k)^(-5)*exp(-B/(W(k)^4));
        S=[S;Sa];
    end
    TitleStr='ITTC-Modified Pierson-Moskowitz Spectrum';
    L1Str = ['Hs =',num2str(Hs),' [m],','  T0 =',num2str(T0),' [s]'];        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
case 4,%ITTC-Modified Pierson-Moskowitz (Hs,T1)
	Hs = Par(1);
    T1 = Par(2); 
    A=173*Hs^2/T1^4;
    B=691/T1^4;
    for k=1:1:length(W),
        Sa=A*W(k)^(-5)*exp(-B/(W(k)^4));
        S=[S;Sa];
    end
    TitleStr='ITTC-Modified Pierson-Moskowitz Spectrum';
    L1Str = ['Hs =',num2str(Hs),' [m],','  T1 =',num2str(T1),' [s]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
case 5,%ITTC-Modified Pierson-Moskowitz (Hs,Tz)
	Hs = Par(1);
    Tz = Par(2); 
    A=123*Hs^2/Tz^4;
    B=495/Tz^4;
    for k=1:1:length(W),
        Sa=A*W(k)^(-5)*exp(-B/(W(k)^4));
        S=[S;Sa];
    end
    TitleStr='ITTC-Modified Pierson-Moskowitz Spectrum';
    L1Str = ['Hs =',num2str(Hs),' [m],','  Tz =',num2str(Tz),' [s]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
case 6,%JONSWAP (Vwind10,Fetch):
	Vw10  = Par(1); 
    fetch = Par(2);
    g=9.81;
    xtilde= g*fetch/(Vw10^2);
    f0=3.5*(g/Vw10)*xtilde^-0.33;
    w0=2*pi*f0;
    alpha=0.076*xtilde^-0.22;
    for k=1:1:length(W),
            if  W(k) < w0,
               sigma=0.07;
            else
               sigma=0.09;
            end
            S1=alpha*g^2*(W(k)^-5)*exp(-(5/4)*(w0/W(k))^4);
            S2=3.3^(exp(-(W(k)-w0)^2/(2*(sigma*w0)^2)));
            S=[S;S1*S2];            
    end
    TitleStr='JONSWAP Spectrum';
    L1Str = ['Vwind @10m ASL =',num2str(Vw10),' [m/s],','  Fetch =',num2str(fetch/1000),' [km]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 7,%JONSWAP (gamma,Hs,w0):
    Hs	  = Par(1);
    w0    = Par(2);
    gamma = Par(3);
    g=9.81;
    alpha=0.2*Hs^2*w0^4/g^2;
    if gamma<1 | gamma > 7,  %DNV recommeded values when gamma is unknown.
		if gamma ~= 0
			% Display warning if gamma is outside validity range and not
			% set to zero
			disp('Warning: gamma value in wave_spectrum function outside validity range, using DNV formula')
		end
		k=2*pi/(w0*sqrt(Hs));
       if k <= 3.6
           gamma = 5;
       elseif k <= 5
		   gamma = exp(5.75-1.15*k);
	   else % k > 5
		   gamma = 1;
	   end
    end
    for k=1:1:length(W),
            if  W(k) < w0,
                 sigma=0.07;
            else
                 sigma=0.09;
            end
            S1=alpha*g^2*(W(k)^-5)*exp(-(5/4)*(w0/W(k))^4);
            S2=gamma^(exp(-(W(k)-w0)^2/(2*(sigma*w0)^2)));
            % DNV conversion factor from <Environmenatal conditions and environmental
            % loads. April 2007, DNV-RP-C205>            
            Conv_factor =  1-0.287*log(gamma); 
            S=[S;Conv_factor*S1*S2];
    end
    TitleStr='JONSWAP Spectrum';
    L1Str = ['gamma =',num2str(gamma),' Hs =',num2str(Hs),' [m],','  w0=',num2str(w0),' [rad/s]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 8,%Torsethaugen (Hs,w0):
    Hs = Par(1);
    w0 = Par(2);
    N=length(W);
    wmax = max(W);
    S = torset_spec(Hs,w0,W'); % See function below
    TitleStr='Torsethaugen Spectrum';
    L1Str = ['Hs =',num2str(Hs),' [m],','  w0=',num2str(w0),' [rad/s]'];
        
otherwise
    disp('Wrong spectrum type identifier, SpecType=1,2,..,8')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Spectrum
if PlotFlag==1,
	plot(W, S,'r','linewidth',1)
	title(TitleStr)
    legend(L1Str)
	xlabel('\omega [rad/s]')
	ylabel('S(\omega) [m^2 s]')
end            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function wave_spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = torset_spec(Hs,wo,omg)
%
% The Torsethaugen spectrum is an empirical two peaked spectrum for swell and 
% developing sea based on experimental data from the North Sea. For small peak
% frequencies, i.e.  0 < wmax <= 0.6 only one peak in the spectrum appears.
% Returns the spectral density function S of the Torsethaugen spectrum for 
% the frequencies:  0 < w < wmax (rad/s).
%
% Ouputs:
%   S     	- vector of power spectral densities (m^2s)
%
% Inputs:
%   Hs    	- significant wave height (m) - mean of the ones third highest waves
%   wo    	- peak frequency (rad/s)
%   omg  	- vector of frequencies at which to calculate S
%
% Ref: K.Torsethaugen (1996): "Model for a Doubly Peaked Wave Spectra"
%      Sintef report no.: STF22 A96204 prepared for Norsk Hydro.
%
% Author:     G. Kleiven, Norsk Hydro 
% Date:       2000-06-15
% Revisions:  2001-07-06,Svein I. Sagatun, Norsk Hydro - minor revisions
%             2001-10-14,Thor I. Fossen - IO compatible with the GNC toolbox
%	          2005-03-12 �yvind Smogeli - Revised to comply with MSS, added output consistency test
%             2007-10-08 �yvind Smogeli - Bug fix for scaling of spectrum magnitude

%---------------------------------------------------------------------------
% source code: Norsk Hydro
%---------------------------------------------------------------------------
Tp   = 2*pi/wo;      % peak period (s)
f2pii = 2*pi;
fwtp = f2pii/Tp;
Nfrq = length(omg);

% Hs must be positive
if Hs > 0
	
	%---------------------------------------------------------------------------
	%  Parameters:
	%---------------------------------------------------------------------------
	af = 6.6;
	ae = 2.0; 
	au = 25;   
	a10 = 0.7;
	a1 = 0.5;
	kg = 35.0;
	kg0 = 3.5;
	kg1 = 1.0;
	r = 0.857;
	k0 = 0.5;
	k00 = 3.2;
	m0 = 4.0;
	b1 = 2.0;
	a20 = 0.6;
	a2 = 0.3;
	a3 = 6.0;
	s0 = 0.08;
	s1 = 3.0;
	b2 = 0.7;
	b3 = 3.0;
	sigma_a2 = 2*0.07^2;
	sigma_b2 = 2*0.09^2;     
	tf = af*Hs^(1./3);
	
	if Tp < tf;     
		%-------------------------------------------------------------------------
		% Predominant wind sea peak:
		%-------------------------------------------------------------------------
		tl = ae*(Hs^0.5);
		eps1 = (tf-Tp)/(tf-tl);
		rpw = (1-a10)*exp(-((eps1/a1)^2))+a10;
		hsw = rpw*Hs;
		hss = sqrt(1.-rpw^2)*Hs;
		tpw = Tp;
		tps = tf+b1;
		sp = ((f2pii/9.81)*hsw/(Tp^2));
		gammaw = kg*(1+kg0*exp(-Hs/kg1))*(sp^r);
		gammas = 1.;
		nw = k0*sqrt(Hs)+k00;
		mw = m0;
		ns = nw;
		ms = mw;
		if ms < 1
			ms = 1;
		end         % Constraint implemented by GKl
		g_argw = (nw-1)/mw;
		g_args = (ns-1)/ms;
		if g_args < 0
			g0s = 1./((1./ms)*gamma(g_args)*((ns/ms)^(-g_args)));
		else
			g0s = 1./((1./ms)*gamma(g_args)/((ns/ms)^(g_args)));
		end
		if g_argw < 0
			g0w = 1./((1./mw)*gamma(g_argw)*((nw/mw)^(-g_argw)));
		else
			g0w = 1./((1./mw)*gamma(g_argw)/((nw/mw)^(g_argw)));
		end
		a1m = 4.1;
		b1m = 2.0*(mw^0.28)-5.3;
		c1m = -1.45*(mw^0.1)+0.96;
		a2m = 2.2/(mw^3.3)+0.57;
		b2m = -0.58*mw^0.37+0.53;
		c2m = -1.04/(mw^1.9)+0.94;
		if c1m < 0 
			f1w = a1m / (nw-b1m)^(-c1m);
		else
			f1w = a1m * (nw-b1m)^c1m;
		end
		if b2m < 0
			f2w = a2m / nw^(-b2m) + c2m;
		else
			f2w = a2m * nw^(b2m) + c2m;
		end
		b1m = 2.0*(ms^0.28)-5.3;
		c1m = -1.45*(ms^0.1)+0.96;
		a2m = 2.2/(ms^3.3)+0.57;
		b2m = -0.58*ms^0.37+0.53;
		c2m = -1.04/(ms^1.9)+0.94;
		if c1m < 0 
			f1s = a1m / (ns-b1m)^(-c1m);
		else
			f1s = a1m * (ns-b1m)^c1m;
		end
		if b2m < 0
			f2s = a2m / ns^(-b2m) + c2m;
		else
			f2s = a2m * ns^(b2m) + c2m;
		end
		
		agammaw = (1+f1w*log(gammaw)^(f2w))/gammaw;
		agammas = (1+f1s*log(gammas)^(f2s))/gammas;
	else
		%----------------------------------------------------------------------------------------------------------------
		% Predominant swell peak:
		%--------------------------------------------------------------------------
		tu = au;
		epsu = (Tp-tf)/(tu-tf);
		rps = (1.-a20)*exp(-(epsu/a2)^2)+a20;
		hss = rps*Hs;
		hsw = sqrt(1.-rps^2)*Hs;
		tps = Tp;
		ns = k0*sqrt(Hs)+k00;
		ms = m0;
		nw = ns;
		mw = m0*(1-b2*exp(-Hs/b3));
		s4 = s0*(1-exp(-Hs/s1));
		g_argw = (nw-1)/mw;
		g_args = (ns-1)/ms;
		if g_args < 0
			g0s = 1./((1./ms)*gamma(g_args)*((ns/ms)^(-g_args)));
		else
			g0s = 1./((1./ms)*gamma(g_args)/((ns/ms)^(g_args)));
		end
		if g_argw < 0
			g0w = 1./((1./mw)*gamma(g_argw)*((nw/mw)^(-g_argw)));
		else
			g0w = 1./((1./mw)*gamma(g_argw)/((nw/mw)^(g_argw)));
		end
		tpw = ((g0w*hsw^2)/(16*s4*(0.4^nw)))^(1./(nw-1.));
		sf = ((f2pii/9.81)*Hs/(tf^2));
		gammaw = 1.;
		gamma_f = kg*(1+kg0*exp(-Hs/kg1))*sf^r;
		gammas = gamma_f*(1.+a3*epsu);
		a1m = 4.1;
		b1m = 2.0*(mw^0.28)-5.3;
		c1m = -1.45*(mw^0.1)+0.96;
		a2m = 2.2/(mw^3.3)+0.57;
		b2m = -0.58*(mw^0.37)+0.53;
		c2m = -1.04/(mw^1.9)+0.94;
		if c1m < 0
			f1w = a1m / (nw-b1m)^(-c1m);
		else
			f1w = a1m * (nw-b1m)^c1m;
		end
		if b2m < 0
			f2w = a2m / nw^(-b2m) + c2m;
		else
			f2w = a2m * nw^(b2m) + c2m;
		end
		b1m = 2.0*(ms^0.28)-5.3;
		c1m = -1.45*(ms^0.1)+0.96;
		a2m = 2.2/(ms^3.3)+0.57;
		b2m = -0.58*(ms^0.37)+0.53;
		c2m = -1.04/(ms^1.9)+0.94;
		if c1m < 0
			f1s = a1m / (ns-b1m)^(-c1m);
		else
			f1s = a1m * (ns-b1m)^c1m;
		end
		if b2m < 0
			f2s = a2m / ns^(-b2m) + c2m;
		else
			f2s = a2m * ns^(b2m) + c2m;
		end
		agammaw = (1+f1w*log(gammaw)^(f2w))/gammaw;
		agammas = (1+f1s*log(gammas)^(f2s))/gammas;
	end
	
	fdenorm_s = (tps*(hss^2))/16;
	fdenorm_w = (tpw*(hsw^2))/16;
	
	%===================================================================
	%  Estimates spectral density for each frequency in array omg:
	%===================================================================
	
	f = omg/f2pii;
	%-------------------------------------------------------------------
	%Wind sea contribution:
	%-------------------------------------------------------------------
	fnw = f*tpw;
	in = max(find(fnw < 1));
	ftest1(1:in) = exp(-(((fnw(1:in)-1).^2)/sigma_a2));
	ftest1(in+1:Nfrq) = exp(-(((fnw(in+1:Nfrq)-1).^2)/sigma_b2));
	gamma_wf = gammaw.^ftest1;
	gamma_ws_1 = fnw.^(-nw);
	gamma_ws_2 = exp(-(nw/mw)*fnw.^(-mw));
	gamma_ws = gamma_ws_1.*gamma_ws_2;
	sw = g0w*agammaw*gamma_ws.*gamma_wf*fdenorm_w/(2*pi);
	%----------------------------------------------------------------------
	%Swell contribution:
	%----------------------------------------------------------------------
	fns = f*tps;
	is = max(find(fns < 1));
	ftest2(1:is) = exp(-(((fns(1:is)-1).^2)/sigma_a2));
	ftest2(is+1:Nfrq) = exp(-(((fns(is+1:Nfrq)-1).^2)/sigma_b2));
	gamma_sf = gammas.^ftest2;
	gamma_ss_1 = fns.^(-ns);
	gamma_ss_2 = exp(-(ns/ms)*fns.^(-ms));
	gamma_ss = gamma_ss_1.*gamma_ss_2;
	ss = g0s*agammas*gamma_ss.*gamma_sf*fdenorm_s/(2*pi);
	%-----------------------------------------------------------------------------
	%Estimates spectral-density Sf(m^2*s)
	%-----------------------------------------------------------------------------
	S = sw + ss;  % frq. in Hz
		
else% If Hs <= 0
	
	S = zeros(size(omg));
end

% Check that the output is real, if not the input is beyond the validity range
if sum(imag(S)) ~= 0

	S = zeros(size(omg));
	
	disp('Torsethaugen spectrum input outside validity range in wave_spectrum, complex output set to zero');
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = S';
