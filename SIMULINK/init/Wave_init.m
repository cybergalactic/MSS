function [Zeta_a, Omega, Phase, Wavenum, Psi] = Wave_init(spectrum_type, hs, omega_peak, psi_mean,...
	gamma, spread, depth, nfreq, ndir, energylim, freq_cutoff, dir_cutoff, rand_freq, rand_dir,...
	rand_seed, plot_spectrum, plot_realization, disp_flag)
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% [Zeta_a, Omega, Phase, Wavenum, Psi] = Wave_init(spectrum_type, hs, omega_peak, psi_mean,...
%	gamma, spread, depth, nfreq, ndir, energylim, freq_cutoff, dir_cutoff, rand_freq, rand_dir,...
%   rand_seed, plot_spectrum, plot_realization, disp_flag)
%
% Project:	"Marine System Simulator: MSS"
%
% Abstract:	Define the harmonic wave components to be used in the simulation by
%			extracting components from a given frequency and direction spectrum.
%
% Inputs:	spectrum_type	- Spectrum type: 1 = ITTC, 2 = JONSWAP, 3 = Torsethaugen Doubly Peaked
%			hs				- Significant wave heigth in sea state [m]
%			omega_peak		- Peak frequency in spectrum, 0 or negative for expectation value based on Hs [rad/s]
%			psi_mean		- Mean direction in spectrum [rad]
%			gamma			- Gamma value for JONSWAP spectrum, if using this, [1,7]. Default is 3.3.
%							  0 for using DNV value based on Hs and omega_peak 
%			spread			- Spreading factor for direction spectrum {1,2,3,4,5} 
%			depth			- Average water depth, for calculation of wave numbers
%			nfreq			- Number of frequency components
%			ndir			- Number of wave directions, 1 for long-crested waves
%			energylim		- Ignore wave components with Sj/S < energylim, [0 1], 0 for using all components.
%							  Set to an integer > 1 for using this many of the most important waves.
%			freq_cutoff		- Cutoff frequency of spectrum = cutoff*omega_mean
%			dir_cutoff		- Direction cutoff of spectrum at sides of grid [rad]
%			rand_freq		- 1 for using random frequnencies in each frequency interval
%			rand_dir		- 1 for using random directions in each direction interval
%			rand_seed		- random number seed, applies to all random numbers (phase, random direction, rand frequency)
%			plot_spectrum	- 1 for plotting of spectrum, 0 for no plotting
%			plot_realization- 1 for plotting of surface relaization, 0 for no plotting
%			disp_flag		- 1 for writing initialization info to screen, 0 for none
%
% Outputs:	Zeta_a			- Vector of harmonic wave amplitudes
%			Omega			- Vector of harmonic wave frequencies
%			Phase			- Vector of harmonic wave phases (random)
%			Wavenum			- Vector of harmonic wave numbers
%			Psi				- Vector of harmonic wave directions
%
% Calls:	GNC toolbox: wavespec.m
%
% Author:	Oyvind Smogeli, September 2004
%
% ________________________________________________________________
%
% MSS GNC is a Matlab toolbox for guidance, navigation and control.
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

if disp_flag ~= 0
	
	disp(' ')
	disp('***********************************************************')
	disp('Simulink Waves block initialization info:')
	
end

if spectrum_type < 4
	
	% Check that input makes sense
	
	% Check that spectrum type is an integer
	spectrum_type = round(spectrum_type);
	
	% Negative hs not allowed
	if hs < 0
		hs = 0;
		if disp_flag ~= 0
			disp(' Input Hs negative, set to zero.')
		end	
	end
	
	% Keep psi_mean within [-pi pi]
	psi_mean = rad2pipi(psi_mean);
	
	% Round the spreading factor to the nearest integer
	spread = round(spread);
	
	% Check that spreading factor is within validity range, correct if not
	if spread < 1 
		spread = 1;
		if disp_flag ~= 0
			disp(sprintf(' Spreading factor set to        : %g', spread))
		end
	elseif spread > 5 
		spread = 5;
		if disp_flag ~= 0
			disp(sprintf(' Spreading factor set to        : %g', spread))
		end
	end
	
	if depth < 0.001
		depth = 0.001;
		if disp_flag ~= 0
			disp(sprintf(' Depth set to 0.001 m (minimum).'))
		end
	end
		
	% at least one frequency
	if nfreq <= 0
		nfreq = 1;
		if disp_flag ~= 0
			disp(' Input zero or negative number of frequencies, set to one.')
		end	
	end
	
	% at least one direction
	if ndir <= 0
		ndir = 1;
		if disp_flag ~= 0
			disp(' Input zero or negative number of directions, set to one.')
		end	
	end
	
	if freq_cutoff < 2
		freq_cutoff = 2;
		if disp_flag ~= 0
			disp(' Frequency cutoff too small, set to 2.')
		end			
	end
	
	if dir_cutoff < 0
		dir_cutoff = 0;
		if disp_flag ~= 0
			disp(' Direction cutoff input negative, set to zero.')
		end			
	elseif dir_cutoff > 3*pi/8
		dir_cutoff = 3*pi/8;
		if disp_flag ~= 0
			disp(' Direction cutoff too large, set to 3*pi/8.')
		end			
	end	
	
	% Calculate expectation value of omega_peak if wanted, based on statistical North Sea data (NORSOK)
	if omega_peak <= 0
		
		omega_peak		= 2*pi/(4.883 + 2.68*hs^0.54);
		
		if disp_flag ~= 0
			
			disp(sprintf(' Calculated peak wave frequency : %0.2f rad/s', omega_peak))
			
		end
		
	end
	
	% Energylim should be a factor between 0 and 1 or an integer > 1
	if energylim < 0
		
		energylim = 0;

		if disp_flag ~= 0
			
			disp(sprintf(' Energylim input negative, set to zero (using all waves in grid)'))
			
		end

		
	elseif energylim > 1
		
		energylim = round(energylim);
		
		% The maximum number of waves is the number of points in the grid
		if energylim > nfreq*ndir
		
			energylim = nfreq*ndir;
			
			if disp_flag ~= 0
				
				disp(sprintf(' energylim > nfreq*ndir, set to nfreq*ndir (using all waves in grid)'))
				
			end

			
		end
		
	end
	
	% Acceleration of gravity
	g = 9.81;

	% Cutoff frequency
	omega_max = freq_cutoff*omega_peak;
	
	% Frequency step
	delta_omega = omega_max/nfreq;
	
	% Frequency vector, starting at delta_omega
	Omega_vec = delta_omega:delta_omega:omega_max;

	% Direction step
	delta_psi = (pi-2*dir_cutoff)/ndir;
		
	% Start direction
	psi_start = psi_mean - pi/2 + delta_psi/2 + dir_cutoff;
	
	% Max direction
	psi_max = psi_mean + pi/2 - delta_psi/2 - dir_cutoff;
	
	% Direction values for loop
	if ndir > 1
		psi_loop = psi_start:delta_psi:psi_max;
	else
		psi_loop = psi_mean;
	end
		
	% Mean wave component energy, total sea state energy = m0 ~= (4*Hs)^2, bounded from below to avoid division by zero
	if hs ~= 0
		mean_energy = (4*hs)^2/(nfreq*ndir);
	else
		mean_energy = 1E-5;
	end
	
	% Wave frequency spectrum, using the GNC toolbox
	switch spectrum_type
		
		case 1 % ITTC modified PM spectrum
			
			% SpecType =3, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T0)
			SpecType = 3;
			Par = [hs, 2*pi/omega_peak];
			
		case 2 % JONSWAP spectrum
			
			% SpecType =7, JONSWAP (p1=Hs,p2=w0,p3=gamma)

			SpecType = 7;
			Par = [hs, omega_peak, gamma];
	
		case 3 % Doubly Peaked spectrum
			

			% SpecType =8, Torsethaugen (p1=Hs,p2=W0) 
			SpecType = 8;
			Par = [hs, omega_peak];

		otherwise
			
			disp('Choose a valid spectrum!!')
		
	end % switch
	
	% Generate the spectral densities, no plotting
	S_vec = wavespec(SpecType,Par,Omega_vec',0);

	% Set first value of Zeta_a vector negative for testing
	Zeta_a(1) = -1;
		
	% Set random generator
	rand('state',rand_seed);
	
	% Evenly distributed random phase
	phase_vector = rand(1,nfreq*ndir)*2*pi;

	% Evenly distributed random directions
	rand_psi_vector = delta_psi*(rand(1,nfreq*ndir) - 1/2);
	
	% Evenly distributed random frequencies
	rand_omega_vector = delta_omega*(rand(1,nfreq*ndir) - 1/2);
	
	% Random number index
	r = 0;
	
	% Harmonic wave component index
	k = 0;
	
	% Frequency index
	m = 0;
	
	% Loop over frequencies
	for omega = Omega_vec
		
		% Step frequency index
		m = m + 1;
		
		% Reset direction index for each frequency
		n = 0;
		
		% Save frequency for plotting
		Omega_plot(m) = omega;
		
		s_omega = S_vec(m);
		
		
		% Loop over directions
		for psi = psi_loop
			
			% Step direction index
			n = n + 1;
			
			% Step random number index
			r = r + 1;
			
			% Calculate direction spectrum value for each direction component
			if ndir > 1
                K = (2^(2*spread-1)*factorial(spread)*factorial(spread-1))/(pi*factorial(2*spread-1));
                s_psi = K*cos(psi-psi_mean)^(2*spread);
            else
                s_psi = 1/delta_psi;
            end
            
			% Keep wave direction within [-pi, pi]
			psi = rad2pipi(psi);
			
			% Save direction for plotting
			Psi_plot(n) = psi;
			
			% Wave component energy
			comp_energy = s_omega*s_psi*delta_omega*delta_psi;

			% Save harmonic wave components if energy is above limit, store all if using the set
			% number of waves approach (energylim > 1, nwaves = energylim)
			
			if comp_energy/mean_energy >= energylim | energylim > 1
				
				% Step component index
				k = k + 1;

				% Calculate amplitude
				Zeta_a(k) = sqrt(2*s_omega*s_psi*delta_omega*delta_psi);

				% Frequency.
				if rand_freq == 1

					% Use randomly chosen frequnencies around the equally spaced ones
					% rand(1) is equally distributed between 0 and 1
					Omega(k) = omega + rand_omega_vector(r);
					
				else

					% Use equally spaced frequencies
					Omega(k) = omega;

				end

				% Wave number
				if depth == inf
					
					% Deepwater dispersion relation
					Wavenum(k) = Omega(k)^2/g;
				
				else
					
					% Find wave number from finite water depth dispersion relation by iteration
					wavenum_0 = Omega(k)^2/g;
					Wavenum(k) = fzero(@(w) w*tanh(w*depth) - Omega(k)^2/g,[wavenum_0,1E10]);
					
				end
				% Direction
				if ndir > 1 & rand_dir == 1

					% Use randomly chosen directions around the equally spaced ones
					% Do not use this for ndir == 1, since this is long-crested waves in one direction
					Psi(k) = psi + rand_psi_vector(r);

				else

					% Use equally spaced directions
					Psi(k) = psi;

				end
			
				% Phase
				Phase(k) = phase_vector(r);		
				
				% Spectrum value
				S_tot(k) = s_omega*s_psi;

				% Spectrum value (for plotting)
				S_plot(m,n) = s_omega*s_psi;

			else % Component discarded due to low energy

				% Spectrum value for plotting set to zero
				S_plot(m,n) = 0;

			end % Energy test

		end % direction loop
		
	end % frequency loop
	
	% Store the total number of waves
	nwaves = k;
		
	% If using the set number of waves approach (energylim < 0, nwaves = - energylim)
	if energylim > 1 
		
		nwaves = energylim;
		
		% Sort the waves according to energy content
		[S_tot,index] = sort(S_tot);
		Zeta_a = Zeta_a(index);
		Psi = Psi(index);
		Omega = Omega(index);
		Wavenum = Wavenum(index);
		Phae = Phase(index);
		
		largewaves = (length(S_tot) - nwaves + 1):(length(S_tot));
			
		% Keep only the nwaves largest components
		S_tot = S_tot(largewaves);
		Zeta_a = Zeta_a(largewaves);
		Psi = Psi(largewaves);
		Omega = Omega(largewaves);
		Wavenum = Wavenum(largewaves);
		Phase = Phase(largewaves);
		
	end
	
	% Return error if energy limit set too high such that no waves are stored
	if Zeta_a(1) < 0
		
		if hs > 0
			
			if disp_flag > 0.5
				
				disp(' ')
				disp(' Wave Initialization Error: Wave component energy limit set too high, no wave components generated')
				disp(' ')
				
			end
			
			% Output zero
			Zeta_a	= 0;
			Omega	= 0;
			Phase	= 0;
			Wavenum	= 0;
			Psi		= 0;
			
		% Only one wave component returned for hs = 0
		else
			
			% Output zero
			Zeta_a	= 0;
			Omega	= 1;
			Phase	= 0;
			Wavenum	= 1/g;
			Psi		= psi_mean;
			
		end
		
		
	else
		
		if disp_flag ~= 0
			
			disp([' Wave components in grid        : ', num2str(nfreq*ndir)])
			disp([' Used wave components           : ', num2str(length(Zeta_a))])
			disp([' Removed wave components        : ', num2str(nfreq*ndir - length(Zeta_a))])
			
		end
		
		% Plot spectrum
		if plot_spectrum > 0
			
			if ndir > 1
				
				[theta, r] = meshgrid(Psi_plot,Omega_plot);
				[X,Y] = pol2cart(theta, r);
				
				figure
				set(gcf,'Position',[230   530   560   420])
				handle_1 = surf(Y,X,S_plot);
				shading interp % use interpolated surface coloring
				set(handle_1, 'FaceAlpha', [0.5]) % set the surface transparent
				xlabel('y')
				ylabel('x')
				zlabel('S [m^2s]')
				title(sprintf('Wave spectrum, mean direction %0.0f degrees',psi_mean*180/pi))
				set(gca,'XTickLabelMode','manual')
				set(gca,'XTickLabel',{' '})
				set(gca,'YTickLabelMode','manual')
				set(gca,'YTickLabel',{' '})
				set(handle_1,'EdgeColor','Black')

				figure
				set(gcf,'Position',[800   530   560   420])
				handle_2 = mesh(Psi_plot*180/pi,Omega_plot,S_plot);
				hold on
				plot3(Psi*180/pi, Omega, S_tot, 'r*')
				set(handle_2,'FaceColor','none') % set no surface color
				xlabel('Dir [deg]')
				ylabel('\omega [rad/s]')
				zlabel('S [m^2s]')
				title(sprintf('Wave spectrum with harmonic components as red stars, mean direction %0.0f degrees',psi_mean*180/pi))
				
			elseif ndir == 1
				
				figure
				plot(Omega_plot, S_plot)
                hold on
                plot(Omega,S_tot, 'r*')
				xlabel('\omega (rad/s)')
				ylabel('S [m^2s]')
				title(sprintf('Long-crested wave spectrum, direction %0.0f degrees, harmonic components as red stars',psi_mean*180/pi))
				
				
			end
			
		end % plot
		
		if plot_realization > 0
			
			% Plot realization of sea surface. Time-consuming!
			
			xdist = 500;
			ydist = 500;
			
			xpoints = 150;
			ypoints = 150;
			
			t = 0;
			
			az = -16;
			el = 80;
			
			%----------------------------------------------------------------
			
			Zeta = zeros(xpoints+1, ypoints+1);
			Z_temp = Zeta;
			
			dx = xdist/xpoints;
			dy = ydist/ypoints;
			
			xvec = 0:dx:xdist;
			yvec = 0:dy:ydist;
				
			hh = waitbar(0,'Running sea state realization');
			
			xc = 0;
								
			for x = xvec
					
				xc = xc + 1;
					
				yc = 0;
					
				for y = yvec

					yc = yc + 1;
		
					Zeta(xc,yc) = sum(Zeta_a.*cos(Omega*t + Phase - Wavenum.*(x*cos(Psi) + y*sin(Psi))));

				end
	
				waitbar(x/length(xvec), hh);
			
			end
			
			close(hh)
			
			% Test if hs is correct
			hs_sim=4*sqrt(mean(mean(Zeta.^2)));
			
			if disp_flag ~= 0
				
				disp(' ')
				disp(' Sea state total energy test:')
				disp(['  hs input : ', num2str(hs), ' m'])
				disp(['  hs output:   ', num2str(hs_sim), ' m'])
				
			end
			
			figure
			set(gcf,'Position',[470   100   560   420])
			mesh(yvec,xvec,Zeta)
			xlabel('y [m]')
			ylabel('x [m]')
			zlabel('z [m]')
			title(sprintf('Sea state realization, %d wave components',nwaves))
			view(az,el)
			
			
		end
		
	end % if exist Zeta_a
	
	
else
	
	if disp_flag ~= 0
		
		disp(' Invalid spectrum type, no waves generated')
		
	end
	
	% Output zero
	Zeta_a	= 0;
	Omega	= 0;
	Phase	= 0;
	Wavenum	= 0;
	Psi		= 0;
	
end % if spectrum_type < 4

if disp_flag ~= 0
	
	disp('***********************************************************')
	
end

